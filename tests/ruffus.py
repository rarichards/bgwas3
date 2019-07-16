from ruffus import *

def main():
    starting_files = ["a.fasta", "b.fasta", "c.fasta"]

    for ff in starting_files:
        open(ff, "w")

    @transform(starting_files,                     
                suffix(".fasta"),                  
                ".sam")                            
    def map_dna_sequence(input_file,
                        output_file):
        ii = open(input_file)
        oo = open(output_file, "w")

    @transform(map_dna_sequence,                   
                suffix(".sam"),                    
                ".bam")                            
    def compress_sam_file(input_file,
                          output_file):
        ii = open(input_file)
        oo = open(output_file, "w")

    @transform(compress_sam_file,                  
                suffix(".bam"),                    
                ".statistics",                     
                "use_linear_model")                
    def summarise_bam_file(input_file,
                           output_file,
                           extra_stats_parameter):
        ii = open(input_file)
        oo = open(output_file, "w")

    pipeline_run()

if __name__=='__main__':
    main()
