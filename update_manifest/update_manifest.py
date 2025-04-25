import os
import yaml
from google.cloud import storage
from datetime import datetime, timezone
from collections import defaultdict
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Generate a YAML manifest of a GCS bucket.")
    parser.add_argument("bucket", type=str, help="Name of the GCS bucket to generate a manifest for")
    parser.add_argument("-o", "--output", type=str, default="manifest.yaml", help="Output YAML file path (default: 'manifest.yaml')")
    parser.add_argument("-m", "--max-files", type=int, default=5, help="Max number of files before a directory is collapsed (default: >= 5)")
    parser.add_argument("-c", "--collapse-all", action="store_true", help="Collapse all directories if max_files is met or exceeded (default: False)")
    parser.add_argument("-d", "--directory-to-collapse", type=str, default=None, help="Collapse all files in indicated directory if max_files is met or exceeded (default: None)")
    return parser.parse_args()


def human_readable_size(size_in_bytes):
    size_in_kb = size_in_bytes / 1024
    if size_in_kb < 1:
        return "{}B".format(size_in_bytes)
    if size_in_kb < 1000:
        return "{}KB".format(round(size_in_kb, 2))
    size_in_mb = size_in_kb / 1024
    if size_in_mb < 1000:
        return "{}MB".format(round(size_in_mb, 2))
    size_in_gb = size_in_mb / 1024
    return "{}GB".format(size_in_gb, 2)

def build_tree(blob_list, bucket, max_files, collapse_all, directory_to_collapse):
    root = {
        'name': bucket,
        'path': f'gs://{bucket}/',
        'type': 'bucket',
        'manifest_generated_date': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M (UTC)'),
        'children': {}
    }

    dir_to_blobs = defaultdict(list)
    dir_sizes = defaultdict(int)

    for blob in blob_list:
        parts = blob.name.strip("/").split('/')
        current_node = root['children']

        for i, part in enumerate(parts):
            is_file = False 
            if (i == len(parts) - 1 and not blob.name.endswith('/')):
              # this is true if it is the last part of the filepath (len(parts) - 1)) and not a directory
              is_file = True
              
            full_path = '/'.join(parts[:i+1])

            if part not in current_node:
                node = {
                    'path': f"gs://{bucket}/{full_path}",
                    'type': 'file' if is_file else 'directory',
                    'last_modified': blob.updated.strftime('%Y-%m-%dT%H:%M (UTC)'),
                }
                if is_file:
                    node['size'] = human_readable_size(blob.size)
                    node['content_type'] = blob.content_type
                else:
                    node['children'] = {}
                current_node[part] = node

            if is_file:
                parent_path = '/'.join(parts[:-1])
                dir_sizes[parent_path] += blob.size
                dir_to_blobs[parent_path].append(blob)
            else:
                current_node = current_node[part]['children']

    def collapse_large_dirs(node, path=""):
        if 'children' not in node:
            return

        for key, child in list(node['children'].items()):
            if child['type'] == 'directory':
                subpath = child['path'].replace(f'gs://{bucket}/', '').strip('/')
                collapse_large_dirs(child, subpath)

        all_children_are_files = all(
            child.get('type') != 'directory'
            for child in node.get('children', {}).values()
        )

        if all_children_are_files and len(node['children']) >= max_files:
            should_collapse = (
                collapse_all or
                (directory_to_collapse and path.startswith(directory_to_collapse))
            )

            if should_collapse:
                blobs = dir_to_blobs.get(path, [])
                if blobs:
                    total_size = sum(blob.size for blob in blobs)
                    last_modified = max(blob.updated for blob in blobs)
                    node['children'] = {
                        'content_summary': {
                            'path': f'gs://{bucket}/{path}/',
                            'type': 'directory_summary',
                            'file_count': len(blobs),
                            'total_size': human_readable_size(total_size),
                            'last_modified': last_modified.strftime('%Y-%m-%dT%H:%M (UTC)'),
                            'note': f"This directory contains more than {max_files} file(s) and has been collapsed into a summary."
                        }
                    }

    def update_directory_sizes(node):
        if node['type'] != 'directory' and node['type'] != 'bucket':
            # Return size in bytes if it's a file
            size_str = node.get('size', '0B')
            if size_str.endswith('KB'):
                return float(size_str[:-2]) * 1024
            elif size_str.endswith('MB'):
                return float(size_str[:-2]) * 1024 * 1024
            elif size_str.endswith('GB'):
                return float(size_str[:-2]) * 1024 * 1024 * 1024
            elif size_str.endswith('B'):
                return float(size_str[:-1])
            return 0

        total = 0
        for child in node.get('children', {}).values():
            total += update_directory_sizes(child)

        node['sum_of_content_size'] = human_readable_size(total)
        return total
    
    update_directory_sizes(root)
    
    if collapse_all or directory_to_collapse != None:
        collapse_large_dirs(root)
  
    return root

def generate_manifest(bucket_name, output_file, max_files, collapse_all, directory_to_collapse):
  
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = list(client.list_blobs(bucket))

    manifest = build_tree(blobs, bucket_name, max_files, collapse_all, directory_to_collapse)

    with open(output_file, 'w') as f:
        yaml.dump(manifest, f, default_flow_style=False, sort_keys=True)

    print(f"Manifest saved to {output_file}")

if __name__ == '__main__':
    args = parse_args()
    generate_manifest(
        bucket_name=args.bucket, 
        output_file=args.output, 
        max_files=args.max_files,
        collapse_all=args.collapse_all,
        directory_to_collapse=args.directory_to_collapse.strip().strip('/')
    )