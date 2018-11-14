require 'net/http'
require 'json'


class Protein

  attr_accessor :gene_id       #Gene id
  attr_accessor :accession     #Array of accessions
  attr_accessor :kegg          #Array with kegg reaction ids
  attr_accessor :kegg_pathway  #Hash kegg reaction id => pathway
  attr_accessor :interactions  #Array with the accessions of the interactors
  attr_accessor :go            #Hash go id => go term
  attr_accessor :intact        #Array of accessions for IntAct
  @@my_proteins = []           #Array with all the proteins
  @@my_intacts = []            #Array with all the proteins with entries on IntAct
  
  def initialize (params = {})
    
    @gene_id = params.fetch(:gene_id, nil)
    #Check if the gene id has the correct format
    unless @gene_id =~/^A[Tt]\d[Gg]\d\d\d\d\d$/
      abort ("The gene ID #{@gene_id} is not valid")
    end
    
    @interactions = []
    @go = {}
    @kegg = []
    @intact = []
    @kegg_pathway = {}
    self.annotate
    
    @@my_proteins << self
    
  end
  
  def annotate
    #Using the gene id annotates the protein (accessions, KEGG reactions id, GO terms and IntAct accessions)
    
    address = URI("http://togows.dbcls.jp/entry/uniprot/#{@gene_id}.json")
    response = Net::HTTP.get_response(address)  #use the Net::HTTP object "get_response" method to call that address
    case response
      when Net::HTTPSuccess then
      #successful retrieval of web page
        data = JSON.parse(response.body)
        
        #Accessions
        @accession = data[0]["accessions"] #Array with all the accessions
        #KEGG
        @kegg =  data[0]["dr"]["KEGG"] #Array of arrays [[kegg_id, -],[],[]...]       
        #GO
        unless data[0]["dr"]["GO"].nil?
          data[0]["dr"]["GO"].each do |annot|
            (id, ont)=annot
            @go[id] = ont if ont.match(/^P:/) #Create a hash with go id => go term (only biological processes)
          end
        end
        #IntAct
        unless data[0]["dr"]["IntAct"].nil?
          data[0]["dr"]["IntAct"].each do |annot|
            @intact << annot[0] #Some proteins have more than 1 accession, so I add them to an array
          end
          @@my_intacts << self #Array with the proteins which have an entry in IntAct
        end

      else
        #unsuccessful retrieval of web page
        raise Exception, "Something went wrong... the call to #{address} failed; type #{response.class}"    
    end
  end 
  
  def annotate_kegg
    #Gets the kegg pathways
    if @kegg
      address = URI("http://togows.dbcls.jp/entry/kegg-reaction/ath:#{@gene_id}.json") #Some of the ids in @kegg didn't correspond to the gene id or A.thaliana, so I only retrieve the pathways for the gene id
      response = Net::HTTP.get_response(address)
      data = JSON.parse(response.body)
      @kegg_pathway = data[0]['pathways']
    end
  end
  
  def Protein.all_intacts
    #Returns an array with all the protein of the list with an entry on IntAct
    return @@my_intacts
  end
  
  def Protein.all_proteins
    #Returns an array with all the objects of the class
    return @@my_proteins
  end
  
  def Protein.find_protein_by_accession(accession)
    #Return a Protein object using its accession
    return @@my_proteins.find { |prot| prot.accession.include? accession }
  end
    
  def Protein.load_from_file(file_name)
    #Load new objects from a file
    if File.file?(file_name)
      File.readlines(file_name).each do |gene_id|
        self.new(:gene_id => gene_id.chomp)
      end
    else
      abort ("Couldn't find the file #{file_name}") #If the file doesn't exist stops the program
    end
  end
  
end