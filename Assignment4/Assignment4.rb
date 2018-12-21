# Assignment 4    Patricia Blanco Gabella
#
# Find orthologs searching Reciprocal Best Hits (RBH) with BLAST
#
# Ortholog genes are defined as genes that are found in different species that evolved from a common
# ancestral gene by speciation. One of the first approaches to detect ortologs is finding RBH. This
# program makes a BLAST using the queries from a genome or proteome against another genome or proteome
# database. After that, it makes another BLAST with the best hits against the first proteome. If they
# are reciprocal, they are probably orthologs.
#
# To make the BLAST I have used the following parameters: e-value <= 10e-10, overlapping >= 50% and
# identity >= 30%. I have read that the e-value depends on the similarity between the organisms and
# on the size of the database. The lower e-value, the more significant is the alignment. Some of the
# papers I read use higher e-value threshold, but I wanted to be restrictive. To minimize the number
# of false positive I also filtered using an overlapping threshold of at least 50% and sequence identity
# had to be larger than 30%. These two values are the most common thresholds I have found in the papers
# I have consulted.
#
# This RBH approach is useful and easy but to be more precise at the time of finding an ortolog pair it
# would be necessary to continue the analysis. One problem of the RBH method is that we could be missing
# orthologs pairs in those cases in which BLAST algorithm returns a paralog as the best hit in at least
# one direction, so the ortholog wont be detected even if it was among the high-scoring hits.It would be
# a good idea to study not only the best hit but the ones with higher scores. Another thing we could try
# is to elaborate a phylogenetic tree to study phylogenetic relationships between the genes.
#
# Sources:
# Huvet, Maxime, and Michael PH Stumpf. "Overlapping genes: a window on gene evolvability." BMC genomics 15.1 (2014): 721.
# Wall, D. P., H. B. Fraser, and A. E. Hirsh. "Detecting putative orthologs." Bioinformatics 19.13 (2003): 1710-1711.
# Kuzniar, Arnold, et al. "The quest for orthologs: finding the corresponding gene across genomes." Trends in Genetics 24.11 (2008): 539-551.


require 'bio'
require 'stringio'

##===================================================================================================
##   FUNCTIONS
##===================================================================================================

def create_db (blast_type,db)
  # Runs makeblastdb to create the database files depending on blast type
  
  case blast_type
  when "blastn", "tblastn"
    %x[makeblastdb -in #{db} -dbtype 'nucl']
    return "n"
  when "blastp", "blastx"
    %x[makeblastdb -in #{db} -dbtype 'prot']
    return "p"
  else
    puts "Blast type not valid. Ignoring..."
    return false
  end
end


def determine_blast (query,db)
  # Determine blast type depending on the kind of sequences in the query file and the database
  
  query_entry = query.next_entry.to_biosequence
  db_entry = db.next_entry.to_biosequence
  
  if query_entry.guess.equal? Bio::Sequence::NA and db_entry.guess.equal? Bio::Sequence::NA
    blast = 'blastn'
  elsif query_entry.guess.equal? Bio::Sequence::AA and db_entry.guess.equal? Bio::Sequence::AA
    blast = 'blastp'
  elsif query_entry.guess.equal? Bio::Sequence::NA and db_entry.guess.equal? Bio::Sequence::AA
    blast = 'blastx'
  elsif query_entry.guess.equal? Bio::Sequence::AA and db_entry.guess.equal? Bio::Sequence::NA
    blast = 'tblastn'
  else
    puts "Couldn't determine blast type"
    blast = nil
  end
  return blast
end


def find_RBH(file1,file2)
  # Compares the results of both blast and looks for RBH
  
  flatfile1 = Bio::FlatFile.auto(file1)
  flatfile2 = Bio::FlatFile.auto(file2)
  blast_1_2 = determine_blast(flatfile1,flatfile2) # BLAST type of the first search
  blast_2_1 = determine_blast(flatfile2,flatfile1) # BLAST type of the second search
  
  return false if blast_1_2 == false or blast_2_1 == false # If the blast type could not be determined, does not make the blast
  
  best_hits = {}   # Hash with the best hits of the first search
  orthologs = {}   # Hash wirh the pairs of orthologs
  
  type_db2 = create_db(blast_1_2,file2) # Create the database for the first blast and save the type (n for nucleotide and p for amino acid)
  flatfile1.rewind
  
  i = 0
  
  flatfile1.each_entry do |query|
    
    puts "Searching... #{i}  #{query.entry_id}"
    best_hit_id = run_BLAST(query,file2,blast_1_2) # Run BLAST for each query
    best_hits[query.entry_id] = best_hit_id if best_hit_id != false
    
    i += 1
  end
  
  type_db1 = create_db(blast_2_1,file1) # Create the database for the secpnd blast and save the type (n for nucleotide and p for amino acid)
  
  flatfile2.rewind
  flatfile2.each_entry do |query|
    if best_hits.values.include? query.entry_id
      puts "Searching... #{query.entry_id}"
      best_hit_id = run_BLAST(query,file1,blast_2_1) # Run BLAST for each query
      
      if best_hits[best_hit_id] == query.entry_id # Check if both results are reciprocal
        orthologs[query.entry_id] = best_hit_id   # Save the ids in the hash if they are RBH
      end
    else # Dont run BLAST with those genes that are not best hits of the other genome
      next
    end
  end
  
  # Delete database files
  %x[rm #{file1}.#{type_db1}hr]
  %x[rm #{file1}.#{type_db1}in]
  %x[rm #{file1}.#{type_db1}sq]  
  %x[rm #{file2}.#{type_db2}hr]
  %x[rm #{file2}.#{type_db2}in]
  %x[rm #{file2}.#{type_db2}sq]
  
  return orthologs
    
end

def run_BLAST(query,db,blast_type)
  # Run BLAST using the desired parameters
  
  evalue = 10e-10
  coverage_threshold = 0.5
  identity_threshold = 0.3
  
  factory = Bio::Blast.local(blast_type, db, "-e #{evalue}")
  
  report = factory.query(query) 
  unless report.hits[0] == nil 
    
    # Calculate coverage of the first hit
    case blast_type
    when "blastx"
      coverage = report.hits[0].overlap.to_f/(report.hits[0].query_len/3).to_f
    when "blastp","blastn","tblastn"
      coverage = report.hits[0].overlap.to_f/report.hits[0].query_len.to_f
    end
    
    identity = report.hits[0].identity.to_f/report.hits[0].overlap.to_f # Calculate identity
    
    # If the coverage and the identity % are above the threshold returns the id of the first hit
    return report.hits[0].definition.split("|")[0].strip if coverage >= coverage_threshold and identity >= identity_threshold 
  end
    
  return false # If there are no hits or they aren't good enough return false
  
end

##===================================================================================================
##   MAIN
##===================================================================================================

usage = "usage: ruby Assignment4.rb file1.fa file2.fa"
abort (usage) if ARGV.length != 2 # Check if both proteomes are present
abort ("Couldn't find the file #{ARGV[0]}") unless File.file?(ARGV[0]) # Check if the file exists
abort ("Couldn't find the file #{ARGV[1]}") unless File.file?(ARGV[1]) # Check if the file exists

file1 = ARGV[0]
file2 = ARGV[1]

orthologs = find_RBH(file1,file2)
abort ("Couldn't find ortholog genes") if orthologs == false

out = File.open("RBH_results.txt", 'w')
orthologs.each {|key, value| out.puts "#{key}\t#{value}" }
