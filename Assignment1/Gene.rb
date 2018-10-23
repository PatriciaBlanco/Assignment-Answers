class Gene

  attr_accessor :id #gene id
  attr_accessor :name #gene name
  attr_accessor :mutant #mutant phenotype
  @@my_genes = [] #array with all gene objects
  @linked = nil #gene object to which it is linked
  
  def initialize (params = {})
    @id = params.fetch(:id, "000")
    #Check if the gene id has the correct format
    unless @id =~/^A[Tt]\d[Gg]\d\d\d\d\d$/
      abort ("The gene ID #{@id} is not valid")
    end
    @name = params.fetch(:name,'Unknown name')
    @mutant = params.fetch(:mutant, 'Unknown phenotype')
    #Add the new object to the array
    @@my_genes << self
  end
  
  def Gene.all_genes
    #Returns an array with all the objects of the class
    return @@my_genes
  end

  def Gene.load_from_file(file_name)
    #Load the new objects from a file
    if File.file?(file_name)
      File.readlines(file_name).drop(1).each do |line|
        id, name, mutant = line.split("\t")
        self.new(:id => id, :name => name, :mutant => mutant)
      end
    else
      abort ("Couldn't find the file #{file_name}") #If the file doesn't exist stops the program
    end
  end

  
  def linked
    return @linked
  end
  
  def linked=(newvalue)
    @linked = newvalue
  end
  
end