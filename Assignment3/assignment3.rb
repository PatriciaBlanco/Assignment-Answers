# Assignment 3    Patricia Blanco Gabella
# For each gene contained in the provided file it creates a Bio::EMBL object. Then, it and looks for CTTCTT repeats
# within its exons and saves them as Bio::Feature objects in the corresponding Bio::Sequence object of the gene.
# Finally it prints this new features in two gff3 files, one with positions referred to the gene, and other referred
# to the chromosome. It also prints a report with the genes without CTTCTT repeats.

require 'bio'
require 'net/http'

#==============Functions==============

def add_repeats (entry, gene_id)
	# Looks for repeats in exons and adds the feature with the position to the Bio::Sequence object
	# Returns the Bio::Sequence with the new repeat features or false if no repeat was found
	
	pattern_f = Bio::Sequence::NA.new("ctt") # Repeat sequence in forward strand
	pattern_c = pattern_f.reverse_complement # Repeat sequence in complement strand
	n_rep = 2 # cttctt = 2 repetitions of the pattern
	len_rep = (pattern_f*n_rep).length # Total length of the repeat
	
	return false unless entry.seq.match(pattern_f*n_rep) or entry.seq.match(pattern_c*n_rep)	# If the sequence doesn't even contain the target sequence return false
	
	bioseq = entry.to_biosequence	# Convert Bio::EMBL object to Bio::Sequence object
	added_pos = [] # Array with the positions of the repeats of the gene that have already been added
	
	entry.ft do |feature|	
		feature_type = feature.feature
		next unless feature_type == "exon" # If the feature is not an exon pass to next feature
		locations_obj = feature.locations # Bio::Locations object that contains Bio::Location objects. Usually it has only one Bio::Location object, but it can contain more
		locations_obj.each do |location|
			# Each Bio::Location object contains the position of the exon
			next if location.xref_id # If the position is referred to a remote entry pass 
			from = location.from - 1		# Start of the exon
			to = location.to - 1				# End of the exon
			exon = entry.seq[from..to]	# Exon sequence
			
			case location.strand
			when 1 # Forward strand
				match = exon.scan(pattern_f) do # Search the motif in the exon 
					if exon[$~.offset(0)[0],len_rep] == pattern_f*n_rep # If next nucleotides are ctt save position
						pos = "#{$~.offset(0)[0]+from+1}..#{$~.offset(0)[0]+from+len_rep}"
						unless added_pos.include? pos
							# To avoid duplicates. Add new feature unless the repeat position has already been added
							bioseq.features << create_feature_repeat(pos, gene_id, location.strand)
							added_pos << pos
						end
					end
				end
				
			when -1 # Complement strand
				match = entry.seq[from..to].scan(pattern_c) do
					if exon[$~.offset(0)[0],len_rep] == pattern_c*n_rep # If next nucleotides are agg save position
						pos = "complement(#{$~.offset(0)[0]+from+1}..#{$~.offset(0)[0]+from+len_rep})"
						unless added_pos.include? pos
							# To avoid duplicates. Add new feature unless the repeat position has already been added
							bioseq.features << create_feature_repeat(pos, gene_id, location.strand)
							added_pos << pos
						end
					end				
				end
			end
		end
	end
	
	if added_pos.empty? 
		return false # If no repeat has been found returns false
	else
		return bioseq
	end
end

def create_EMBL(gene_id)
	#Returns Bio::EMBL object using gene_id
	
	response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id}") # Retrieve entry for that gene id
	if response == false 
		puts "Error retrieving gene #{gene_id}. Ignoring..."
		return false #Unsuccessful retrieval
	end

	record = response.body
	return Bio::EMBL.new(record) # Create Bio::EMBL object
end

def fetch(uri_str)
	# Access the provided address and return the results
	
	address = URI(uri_str)
  response = Net::HTTP.get_response(address) # Uses the Net::HTTP object "get_response" method to call the provided address

  case response 
    when Net::HTTPSuccess then 
      # Successful retrieval of web page
      return response
    else
			# Unsuccessful retrieval of web page
      raise Exception, "Something went wrong... the call to #{address} failed; type #{response.class}"
      return false
  end
end

def create_feature_repeat(pos, gene_id, strand)
	# Add new cttctt repeat feature to a Bio::Sequence object
	
	f = Bio::Feature.new('cttctt_repeat', pos)
	f.append(Bio::Feature::Qualifier.new('motif', 'cttctt'))
	if strand == 1
		f.append(Bio::Feature::Qualifier.new('strand', '+'))
	elsif strand == -1
		f.append(Bio::Feature::Qualifier.new('strand', '-'))
	end
	f.append(Bio::Feature::Qualifier.new('gene', gene_id))
	return f
end

def write_gff3_chr (targets)
	# Write a gff3 file with chromosome coordinates
	
	chr_gff3 = File.open("chr_results.gff3", 'w')
	chr_gff3.puts "##gff-version 3\n"
	
	targets.each do |gene_id,bioseq|
		chr_info = bioseq.primary_accession.split(':') # e.g chromosome:TAIR10:4:1229031:1230164:1
		n_id = 1
		bioseq.features.each do |feature|
			featuretype = feature.feature
			next unless featuretype == "cttctt_repeat"
			loc = feature.locations.first
			strand = feature.assoc["strand"]
			my_id = "CTTCTT_repeat_#{gene_id}_#{n_id}"
			
			#seqid	source	start	end	score	strand	phase	attributes
			chr_gff3.puts "Chr#{chr_info[2]}\tAssignment3\trepeat_region\t#{loc.from+chr_info[3].to_i-1}\t#{loc.to+chr_info[3].to_i-1}\t.\t#{strand}\t.\tID=#{my_id};Name=repeat_cttctt"
			n_id += 1
		end			
	end
end

def write_gff3_gene (targets)
	# Write a gff3 file
	
	gff3 = File.open("results.gff3", 'w')
	gff3.puts "##gff-version 3\n"
	
	targets.each do |gene_id,bioseq|
		n_id = 1
		bioseq.features.each do |feature|
			featuretype = feature.feature
			next unless featuretype == "cttctt_repeat"
			loc = feature.locations.first
			strand = feature.assoc["strand"]
			my_id = "CTTCTT_repeat_#{gene_id}_#{n_id}"
			
			#seqid	source	start	end	score	strand	phase	attributes
			gff3.puts "#{gene_id}\tAssignment3\trepeat_region\t#{loc.from}\t#{loc.to}\t.\t#{strand}\t.\tID=#{my_id};Name=repeat_cttctt"
			n_id += 1
		end			
	end
end

#===============================================================================================================

usage = "usage: ruby assignment3.rb ArabidopsisSubNetwork_GeneList.txt"
abort (usage) if ARGV.length != 1 # Check if the gene list file is present
abort ("Couldn't find the file #{ARGV[0]}") unless File.file?(ARGV[0]) # Check if the file exists

print "\nSearching repeats..."
genes = {}			# Hash with all Bio::EMBL objects for the genes in the list (gene_id => Bio::EMBL)
targets = {}		# Hash with the Bio::Sequence objects with the target in exons (gene_id => Bio::Sequence)

File.readlines(ARGV[0]).each do |gene_id|
	gene_id = gene_id.chomp
	entry = create_EMBL(gene_id) # Create a new Bio::EMBL object
	if entry.is_a? Bio::EMBL 
		genes[gene_id] = entry
	else 
		next # Unsuccessful retrieval
	end
	
	bioseq = add_repeats(entry,gene_id) # Bio::Sequence object with the new repeat features or false if no repeat was found
	targets [gene_id] = bioseq unless bioseq == false
	
end
puts "Done\n\n"

#Write gff3 file with the repeats
write_gff3_gene(targets)
write_gff3_chr(targets)

#Write report with genes without repeats
no_rep = File.open("no_repeat.txt", 'w')
genes.keys.each do |gene|
	unless targets.keys.include? gene
		no_rep.puts gene
	end
end

puts "Found repeats have been saved in repeats.gff3 (gene coordinates) and chr_repeats.gff3 (chromosome coordinates)"
puts "Genes without repeats in exons have been saved in no_repeats.txt"
