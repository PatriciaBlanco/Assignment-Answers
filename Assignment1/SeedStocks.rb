require './Gene.rb'

class SeedStocks

  attr_accessor :seed_stock #Seed stock name
  attr_accessor :id #Mutated gene id
  attr_accessor :last_planted #Last planted date
  attr_accessor :storage #Storage
  attr_accessor :grams #Grams remaining
  attr_accessor :gene #Gene object with the same id
  @@my_stocks = [] #Array with all the objects of the class
    
  def initialize (params = {})
    @seed_stock = params.fetch(:seed_stock,'Unknown seed stock')
    @id = params.fetch(:id,"000")
    @last_planted = params.fetch(:last_planted, 'Unknown date')
    @storage = params.fetch(:storage,'Unknown storage')
    @grams = params.fetch(:grams,"000").to_f
    
    #Assign the Gene object with the same id to the variable @gene.
    #If there is no Gene object with that id the code stops
    @gene = Gene.all_genes.find { |gen| gen.id == @id }
    
    if @gene == nil
      abort ("Gene #{@id} associated to seed stock #{@seed_stock} is not registered")
    end
    
    @@my_stocks << self
    
  end

  def plant(some_grams)
    #Substrats some_grams to the grams remaining and updates the date.
    #If the amount of seed is reduced to zero prints a warning message
		if @grams <= 7
      @grams = 0
      puts "WARNING: we have run out of Seed Stock #{@seed_stock}"
    else
      @grams = @grams - some_grams
    end
    @last_planted = Time.now.strftime("%d/%m/%Y")
  end
  
  def SeedStocks.load_from_file(file_name)
    #Load the new objects from a file
    if File.file?(file_name)
      File.readlines(file_name).drop(1).each do |line|
        seed_stock, id, last_planted, storage, grams = line.split("\t")
        self.new(:seed_stock => seed_stock, :id => id, :last_planted => last_planted, :storage => storage, :grams => grams)
      end
    else
      abort ("Couldn't find the file #{file_name}") #If the file doesn't exist stops the program
    end
  end
    
  def SeedStocks.write_database(newfile)
    #Writes a new file with all the updated objects
    file = File.open(newfile, 'w')
    file.write("Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n")
    
    @@my_stocks.each do |stock|
      file.write("#{stock.seed_stock}\t#{stock.id}\t#{stock.last_planted}\t#{stock.storage}\t#{stock.grams.to_i}\n")
    end
    
    file.close
    
  end
  
  def SeedStocks.all_stocks
    #Returns an array with all of the objects of the class
    return @@my_stocks
  end
  
  def SeedStocks.get_seed_stock(seed_stock_name)
    #Returns the object with the provided name
    return @@my_stocks.find { |stock| stock.seed_stock == seed_stock_name }
  end
  
end