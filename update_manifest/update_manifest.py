import yaml
from google.cloud import storage
from datetime import datetime, timezone
from collections import defaultdict
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s  [%(levelname)s] %(message)s')

def parse_args():
    """Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments including bucket name, output file name, 
                            maximum files before collapsing, collapse flag, a 
                            specific directory to collapse, and a previous manifest flag.
    """ 
    parser = argparse.ArgumentParser(description="Generate a YAML manifest of a GCS bucket.")
    parser.add_argument("bucket", type=str, help="Name of the GCS bucket to generate a manifest for")
    parser.add_argument("-o", "--output", type=str, default="manifest.yaml", help="Output YAML file path (default: 'manifest.yaml')")
    parser.add_argument("-m", "--max-files", type=int, default=5, help="Max number of files before a directory is collapsed (default: >= 5)")
    parser.add_argument("-c", "--collapse-all", action="store_true", help="Collapse all directories if max_files is met or exceeded (default: False)")
    parser.add_argument("-d", "--directory-to-collapse", nargs="*", default=[], help="Collapse all files in indicated directories (space-delimited) if max_files is met or exceeded (default: None)")
    parser.add_argument("-p", "--previous-manifest", type=str, default=None, help="Path to a previous YAML manifest file for change detection (default: None)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging output (default: False)")
    return parser.parse_args()

def human_readable_size(size_in_bytes):
    """Converts a size in bytes to human-readable string format (B, KB, MB, GB).
    
    Args:
        size_in_bytes (int): Size in bytes.

    Returns:
        int: Human-readable size string
    """
    size_in_kb = size_in_bytes / 1024
    if size_in_kb < 1:
        return "{}B".format(size_in_bytes)
    if size_in_kb < 1000:
        return "{}KB".format(round(size_in_kb, 2))
    size_in_mb = size_in_kb / 1024
    if size_in_mb < 1000:
        return "{}MB".format(round(size_in_mb, 2))
    size_in_gb = size_in_mb / 1024
    return "{}GB".format(round(size_in_gb, 2))

def convert_size_to_bytes(size_str):
    """Convert a size string (B, KB, MB, GB) into bytes.
    
    Args:
        size_str (str): Size string (e.g., '10KB', '5MB', '2GB').

    Returns:
        int: The size in bytes.
    """
    if size_str.endswith('KB'):
        return int(float(size_str[:-2]) * 1024)
    elif size_str.endswith('MB'):
        return int(float(size_str[:-2]) * 1024 * 1024)
    elif size_str.endswith('GB'):
        return int(float(size_str[:-2]) * 1024 * 1024 * 1024)
    elif size_str.endswith('B'):
        return int(float(size_str[:-1]))
    return 0 

def calculate_size_difference(old_size, new_size):
    """Calculate the size difference between the old and new size.
    
    Args:
        old_size (str): The previous size (as a string with units like 'KB', 'MB', 'GB').
        new_size (str): The new size (as a string with units like 'KB', 'MB', 'GB').

    Returns:
        str: A string representing the size difference with a "+" or "-" sign, or "null" if no change.
    """
    if old_size is None or new_size is None:
        return "null"

    old_size_bytes = convert_size_to_bytes(old_size)
    new_size_bytes = convert_size_to_bytes(new_size)

    size_diff = new_size_bytes - old_size_bytes

    if size_diff > 0:
        return f"+{human_readable_size(size_diff)}"
    elif size_diff < 0:
        return f"{human_readable_size(size_diff)}"
    else:
        return 0

def build_tree(blob_list, bucket, max_files, collapse_all, directory_to_collapse):
    """Build a nested dictionary representing the structure of a GCS bucket.

    Args:
        blob_list (List[storage.Blob]): List of blobs in the GCS bucket.
        bucket (str): Name of the GCS bucket.
        max_files (int): Maximum number of files before collapsing a directory.
        collapse_all (bool): Whether to collapse all large directories.
        directory_to_collapse (str): Specific directory to collapse

    Returns:
        dict: A nested dictionary representing the bucket manifest tree.
    """
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
              # this is true if it is the last part of the filepath (len(parts) - 1)) and not a directory (ends with '/')
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

    def collapse_large_dirs(node, path="", collapse_dir=None):
        """Recursively collapse large directories into summarized entires if they exceed max_files.

        Args:
            node (dict): Current node in the tree.
            path (str): Path of the current node relative to the bucket root.
        """
        if 'children' not in node:
            return

        for key, child in list(node['children'].items()):
            if child['type'] == 'directory':
                subpath = child['path'].replace(f'gs://{bucket}/', '').strip('/')
                collapse_large_dirs(child, subpath, collapse_dir)

        all_children_are_files = all(
            child.get('type') != 'directory'
            for child in node.get('children', {}).values()
        )

        if all_children_are_files and len(node['children']) >= max_files:
            should_collapse = (
                collapse_all or
                (collapse_dir and collapse_dir in path)
            )

            logging.debug("Checking collapse (%s) for path: %s | Should collapse: %s", collapse_dir, path, should_collapse)
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
        """Recursively calculate and annotate the total content size of each directory.

        Args:
            node (dict): Current node in the tree.

        Returns:
            float: Total size of the directory in human-readable format.
        """
        if node['type'] != 'directory' and node['type'] != 'bucket':
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
    
    if collapse_all or len(directory_to_collapse) > 0:
        if len(directory_to_collapse) > 0:
            for directory in directory_to_collapse:
                collapse_large_dirs(root, collapse_dir=directory)
        else:
            collapse_large_dirs(root)
  
    return root

def flatten_manifest(manifest):
    """Flatten the manifest into a dictionary mapping paths to their metadata.

    Args:
        manifest (dict): Manifest tree.

    Returns:
        dict: Flattened path-to-node dictionary.
    """
    flat = {}

    def recurse(node):
        if node.get('type') in ('file', 'directory_summary'):
            flat[node['path']] = node
        elif 'children' in node:
            for child in node['children'].values():
                recurse(child)

    recurse(manifest)
    return flat

def compare_tree_manifests(old_manifest, new_manifest):
    """Compare two flattened manifests to identify added, deleted, and updated files.

    Args:
        old_manifest (dict): Old manifest tree.
        new_manifest (dict): New manifest tree.

    Returns:
        dict: dictionary containing added, deleted, and updated files with their metadata.
    """
    logging.info("Starting comparison of manifests...")
    
    old_flat = flatten_manifest(old_manifest)
    new_flat = flatten_manifest(new_manifest)

    added_files = {}
    deleted_files = {}
    updated_files = {}

    for new_path, new_node in new_flat.items():
        old_node = old_flat.get(new_path)
        new_size = new_node.get('size') or new_node.get('total_size')
        
        if old_node is None:
            added_files[new_path] = {
                "previous_last_modified": None,
                "previous_size": None,
                "size_difference": calculate_size_difference('0B', new_size)
            }
            logging.debug("Added new file/dir: %s", new_path)
            
        else:
            old_size = old_node.get('size') or old_node.get('total_size')
            
            old_last_modified = old_node.get('last_modified')
            new_last_modified = new_node.get('last_modified')

            if old_size != new_size or old_last_modified != new_last_modified:                
                updated_files[new_path] = {
                    "previous_last_modified": old_last_modified,
                    "previous_size": old_size,
                    "size_difference": calculate_size_difference(old_size, new_size)
                }
                
                if 'children' not in new_node:  # if it's a collapsed directory there are no child files
                    logging.debug("Updated directory summary detected: %s (collapsed directory)", new_path)
                    logging.info("Collapsed directory summary updated: %s | Size difference: %s", new_path, calculate_size_difference(old_size, new_size))
                else:
                    logging.debug("Updated file detected: %s", new_path)
    
    for old_path, old_node in old_flat.items():
        if old_path not in new_flat:
            deleted_files[old_path] = {
                "previous_last_modified": old_node.get('last_modified'),
                "previous_size": old_node.get('size') or old_node.get('total_size'),
                "size_difference": None
            }
            logging.debug("Deleted file/dir: %s", old_path)

    changes_summary = {
        "added_files": added_files,
        "deleted_files": deleted_files,
        "updated_files": updated_files,
        "total_number_of_changed_files": len(added_files) + len(deleted_files) + len(updated_files)
    }

    return changes_summary

def generate_manifest(bucket_name, output_file, max_files, collapse_all, directory_to_collapse, previous_manifest_path):
    """Generate and save a YAML manifest of a GCS bucket

    Args:
        bucket_name (str): Name of the GCS bucket.
        output_file (str): Path to save the generated YAML manifest.
        max_files (int): Maximum number of files before collapsing a directory.
        collapse_all (bool): Whether to collapse all directories over the threshold.
        directory_to_collapse (str): Specific directory path to collapse.
        previous_manifest_path (str): Path to a previous manifest for change detection.
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = list(client.list_blobs(bucket))
    if not blobs:
        logging.warning("No blobs found in the bucket: %s", bucket_name)
        return

    manifest = build_tree(blobs, bucket_name, max_files, collapse_all, directory_to_collapse)

    manifest['total_number_of_files'] = len(blobs)
    
    changes = []
    if previous_manifest_path:
        logging.info("Previous manifest provided. Attempting to load and compare: %s", previous_manifest_path)
        try:
            with open(previous_manifest_path, 'r') as f:
                previous_manifest = yaml.safe_load(f)
            changes = compare_tree_manifests(previous_manifest, manifest)
            logging.info("Detected %d changed files.", len(changes))
        except Exception as e:
            logging.error("Failed to load or compare previous manifest: %s", str(e))

    manifest['changed_files'] = changes or {}

    with open(output_file, 'w') as f:
        yaml.dump(manifest, f, default_flow_style=False, sort_keys=False)

    logging.info("Manifest saved to %s", output_file)

if __name__ == '__main__':
    args = parse_args()
    collapse_dirs = set(args.directory_to_collapse)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)
    generate_manifest(
        bucket_name=args.bucket, 
        output_file=args.output, 
        max_files=args.max_files,
        collapse_all=args.collapse_all,
        directory_to_collapse=set(args.directory_to_collapse),
        previous_manifest_path=args.previous_manifest
    )