with open("formatted-nymain.tsv") as infile:
     for line in infile:
             line = line.strip().split('\t')
             for pos in range(int(line[1]), int(line[2])):
                     print(line[0] +'\t'+str(pos) + '\t' + line[3])
 
