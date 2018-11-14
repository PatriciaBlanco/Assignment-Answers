require './Protein.rb'
require './InteractionNetwork.rb'


usage = "usage: ruby assignment2.rb ArabidopsisSubNetwork_GeneList.txt"


if ARGV.length != 1 #Check if the gene list file is present
	puts usage
	exit
else
	puts "Loading genes from file...\n\n"
	Protein.load_from_file(ARGV[0])
end


puts "Creating networks...\n\n"	
Protein.all_intacts.each do |protein|
	#Create the interaction network of the protein objects with an entry on IntAct
  InteractionNetwork.create_network(protein)
end

InteractionNetwork.all_networks.each do |my_network|
	#Annotate the kegg pathway of the proteins in the list that interact in the network
	my_network.protein_list.each {|protein| protein.annotate_kegg }
end

#Write results in a file
file = File.open("Report.txt", 'w')
i = 1
InteractionNetwork.all_networks.each do |my_network|
	file.write("###Network #{i}: #{my_network.protein_list.length} genes from the list are interacting\n")
	file.write("\n#Gene ID\tAccessions\tKEGG pathways\tGO terms\tInteractors")
	my_network.protein_list.each do |int_protein|
		file.write("\n#{int_protein.gene_id}\t#{int_protein.accession}\t#{int_protein.kegg_pathway}\t#{int_protein.go}\t#{int_protein.interactions}")
	end
	#file.write("\n#Network accessions:\n#{my_network.network}\n")
	file.write("\n\n#Network interactions:\n#{my_network.interactions}\n\n\n")
	
	#puts "KEGG: #{my_network.get_kegg_pathway}"
	#puts "GO: #{my_network.get_go_terms}"
	i += 1
end
file.close

puts "FINAL REPORT: #{InteractionNetwork.all_networks.length} interaction networks have been found.\n"
puts "Final results have been saved in Report.txt"
