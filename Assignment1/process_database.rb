require './Gene.rb'
require './SeedStocks.rb'
require './HybridCross.rb'

usage = "usage: ruby process_database.rb  gene_information.tsv  seed_stock_data.tsv  cross_data.tsv  new_stock_file.tsv"

#Check if all the file names are present
if ARGV.length != 4
	puts usage
	exit
else
	gene_information = ARGV[0]
	seed_stock_data = ARGV[1]
	cross_data = ARGV[2]
	new_stock_file = ARGV[3]
end

#Create the objects for each class
Gene.load_from_file(gene_information)
SeedStocks.load_from_file(seed_stock_data)
HybridCross.load_from_file(cross_data)

#Task 1) for each stock 7 grams are planted
SeedStocks.all_stocks.each do |seed|
	seed.plant(7)
end
SeedStocks.write_database(new_stock_file) #Save the updated data in a new file

#Task 2) Calculate the chisquare value for each hybrid cross
HybridCross.all_crosses.each do |cross|
	cross.chi2_test
end

#Bonus: Function to access individual SeedStock objects based on their ID
#puts SeedStocks.get_seed_stock("A334").inspect

#Final report. If the gene is linked to other gene prints the message
puts "\n\nFinal report\n\n"
Gene.all_genes.each do |gen|
	puts "#{gen.name} is linked to #{gen.linked.name}" unless gen.linked == nil
end
