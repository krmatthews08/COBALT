from collections import defaultdict
import sys

def process_blast_output(blast_file):
    # Store hits per query
    query_hits = defaultdict(list)
    
    # Read BLAST results
    with open(blast_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id = fields[0]
            bit_score = float(fields[11])
            query_hits[query_id].append((bit_score, line))
    
    # Print top 5 hits per query
    for query in query_hits:
        sorted_hits = sorted(query_hits[query], reverse=True)
        for score, line in sorted_hits[:1]:
            print(line.strip())

if __name__ == "__main__":
    process_blast_output(sys.argv[1])
