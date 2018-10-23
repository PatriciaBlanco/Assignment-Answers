require './SeedStocks.rb'

class HybridCross

  attr_accessor :p1 #Parent1
  attr_accessor :p2 #Parent2
  attr_accessor :f2_wild #F2 wild type phenotype
  attr_accessor :f2_p1 #F2 parent 1 phenotype
  attr_accessor :f2_p2 #F2 parent 2 phenotype
  attr_accessor :f2_p1p2 #F2 parent 1 and parent 2 phenotype
  attr_accessor :seed1 #Parent 1 seedstock object
  attr_accessor :seed2 #Parent 2 seedstock object
  attr_accessor :chi #Chisquare value for the cross
  @@my_crosses = [] #Array with all the objects in the class
  
  def initialize (params = {})
    @p1 = params.fetch(:p1,'Unknown P2')
    @p2 = params.fetch(:p2,'Unknown P1')
    @f2_wild = params.fetch(:f2_wild, "0").to_f
    @f2_p1 = params.fetch(:f2_p1,"0").to_f
    @f2_p2 = params.fetch(:f2_p2,"0").to_f
    @f2_p1p2 = params.fetch(:f2_p1p2,"0").to_f
    
    #Assign the SeedStocks object with the same name as parent to the variable @seed.
    #If there is no SeedStocks object with that name the code stops
    @seed1 = SeedStocks.all_stocks.find { |seed| seed.seed_stock == @p1 }
    if @seed1 == nil
      abort ("The seed stock #{@p1} is not registered")
    end
    
    @seed2 = SeedStocks.all_stocks.find { |seed| seed.seed_stock == @p2 }
    if @seed2 == nil
      abort ("The seed stock #{@p2} is not registered")
    end
    
    @chi = nil
    
    @@my_crosses << self #Add the new object to the array
    
  end
  
  def HybridCross.all_crosses
    #Returns an array with all the objects of the class
    return @@my_crosses
  end
  
  def HybridCross.load_from_file(file_name)
    #Load the new objects from a file
    if File.file?(file_name)
      File.readlines(file_name).drop(1).each do |line|
        p1, p2, f2_wild, f2_p1, f2_p2, f2_p1p2 = line.split("\t")
        self.new(:p1 => p1, :p2 => p2, :f2_wild => f2_wild, :f2_p1 => f2_p1, :f2_p2 => f2_p2, :f2_p1p2 => f2_p1p2)
      end
    else
      abort ("Couldn't find the file #{file_name}") #If the file doesn't exist stops the program
    end
  end
   
  def chi2_test
    #Calculates the value of chi square for the cross
    
   total = @f2_wild + @f2_p1 + @f2_p2 + @f2_p1p2
   #chisquare = sum((Observed-Expected)²/Expected)
   @chi = (@f2_wild - total*9/16)**2/(total*9/16) + (@f2_p1 - total*3/16)**2/(total*3/16) + (@f2_p2 - total*3/16)**2/(total*3/16) + (@f2_p1p2 - total*1/16)**2/(total*1/16)
   if @chi >= 3.84 #chi distribution value to be considered statistically significant (p<0.05)
     puts "Recording: #{@seed1.gene.name} is genetically linked to #{@seed2.gene.name} with chisquare score #{@chi}"
     @seed1.gene.linked = @seed2.gene
     @seed2.gene.linked = @seed1.gene
   #else
     #puts "No linked with chisquare score #{@chi}"
   end
  end
end