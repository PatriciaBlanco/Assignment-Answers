require './Protein.rb'

class InteractionNetwork
  attr_accessor :network       #Array with the accession of all the proteins in the network
  attr_accessor :origin        #Initial protein object of the network
  attr_accessor :protein_list  #Array with the protein objects that interact in the network
  attr_accessor :interactions  #Hash with the conexions between proteins
  @@my_networks = {}           #Hash with all the filtered networks [interacting proteins of the list] => InteractionNetwork object

  
  def initialize (params = {})
    @network = params.fetch(:network, []) 
    @origin = params.fetch(:origin, nil) #First level of interaction
    @protein_list = params.fetch(:protein_list,[])
    @interactions = params.fetch(:interactions,{})
  end

  def add_interactions(id)
    #Using the accession retrieves the interactors
    arath_taxid = 'taxid:3702(arath)|taxid:3702("Arabidopsis thaliana (Mouse-ear cress)")' #Taxid for A.thaliana
    valid_meth = ['psi-mi:"MI:0006"(anti bait coimmunoprecipitation)',
                  'psi-mi:"MI:0007"(anti tag coimmunoprecipitation)',
                  'psi-mi:"MI:0047"(far western blotting)',
                  'psi-mi:"MI:0055"(fluorescent resonance energy transfer)',
                  'psi-mi:"MI:0065"(isothermal titration calorimetry)',
                  'psi-mi:"MI:0071"(molecular sieving)',
                  'psi-mi:"MI:0084"(phage display)',
                  'psi-mi:"MI:0096"(pull down)',
                  'psi-mi:"MI:0112"(ubiquitin reconstruction)',
                  'psi-mi:"MI:0402"(chromatin immunoprecipitation assay)',
                  'psi-mi:"MI:0663"(confocal microscopy)',
                  'psi-mi:"MI:0676"(tandem affinity purification)',
                  'psi-mi:"MI:1356"(validated two hybrid)'] #Methods I have considered valid

    address = URI("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{id}?format=tab25") #URL to retrieve the interactors in tab25 format
    response = Net::HTTP.get_response(address)  # use the Net::HTTP object "get_response" method to call that address
    case response
      when Net::HTTPSuccess then
      #Successful retrieval of web page
        lines = response.body.split(/\n/)
        interactors = []
        lines.each do |line|
          fields = line.split("\t")
          #score = $1.to_f if fields[14] =~ /intact-miscore:([^"]+)/ #IntAct database regards data with a score of >0.6 as high-confidence and 0.45–0.6 as medium confidence
          if fields[9] == arath_taxid and fields[10] == arath_taxid and valid_meth.include? fields[6] #and score >= 0.45
            #Filter IntAct results depending on the organism, detection methods (and score if required)
            
            id_1 = $1 if fields[0]=~ /uniprotkb:([^"]+)/
            id_2 = $1 if fields[1]=~ /uniprotkb:([^"]+)/
            
            prot = Protein.find_protein_by_accession(id_1)
            unless prot.nil?
              #if the protein belongs to the list, add the interactor to the protein object
              prot.interactions << id_2 unless id_1 == id_2 or prot.interactions.include? id_2 #Don't include the interactor if the protein interacts with itself
              @protein_list << prot unless @protein_list.include? prot 
            end
            
            prot = Protein.find_protein_by_accession(id_2)
            unless prot.nil?
              #if the protein belongs to the list, add the interactor to the protein object
              prot.interactions << id_1 unless id_2 == id_1 or prot.interactions.include? id_1 #Don't include the interactor if the protein interacts with itself
              @protein_list << prot unless @protein_list.include? prot
            end
            
            interactors << id_1 unless id_1.nil? or id_1 == id
            interactors << id_2 unless id_2.nil? or id_2 == id
            
            unless @network.include? id_1 or @origin.accession.include? id_1 or id_1.nil? 
              @network << id_1 #Add new interactor to the network unless it's already in the list or it's another accession for the original protein
            end
            
            unless @network.include? id_2 or @origin.accession.include? id_2 or id_2.nil?
              @network << id_2 #Add new interactor to the network unless it's already in the list or it's another accession for the original protein
            end
          end
        end
        
        @interactions[id] = interactors unless interactors.nil?

      else
        #Unsuccessful retrieval of web page
        puts "No interactions"
        raise Exception, "Something went wrong... the call to #{address} failed; type #{response.class}"
    end
  end
  
  def add_filtered_network
    #Filters new networks
    
    if @@my_networks.keys.include? @protein_list
      #If there's another network with the same interacting proteins, update the network if the new one is smaller (less nodes)
      @@my_networks[@protein_list] = self if @network.length < @@my_networks[@protein_list].network.length
      return
    else
      #The exact key does not match
      @@my_networks.keys.each do |key|
        if (key & @protein_list) == @protein_list #Check if the protein list is included in a bigger network
          return
        elsif (key & @protein_list) == key #Check if the protein list includes a smaller network
          if @protein_list.length > key.length
            @@my_networks.delete(key) #Delete the keys included in @protein_list
          end
        end
      end
      #Add new networks to the hash if the network includes at least 2 proteins from the list
      @@my_networks[@protein_list] = self unless @protein_list.length <= 1 or @@my_networks.keys.include? @protein_list
    end    
  end 
  
  def get_go_terms
    #Gets the go terms of the registered proteins in the network
    go_terms =  {}
    @protein_list.each {|protein| go_terms = go_terms.merge(protein.go)}
    return go_terms
  end
  
  def get_kegg_pathway
    #Gets the kegg_pathway of the registered proteins in the network
    kegg_pathway = {}
    @protein_list.each {|protein| kegg_pathway = kegg_pathway.merge(protein.kegg_pathway)}
    return kegg_pathway
  end

  def InteractionNetwork.all_networks
    #Returns an array of the filtered InteractionNetwork objects
    return @@my_networks.values
  end

  def InteractionNetwork.create_network(origin_protein)
    #Creates a new InteractionNetwork object for a Protein object
    
    #Initialize a new InteractionNetwork object. The first element of the network is the first accession of the origin protein
    new_network = self.new(:network => [origin_protein.intact[0]], :origin => origin_protein)
    
    #2nd level network
    origin_protein.intact.each do |origin_accession| #@intact contains the accession(s) for the protein in IntAct
      new_network.add_interactions(origin_accession) #Obtain interactors of the origin protein using the accession
    end
    
    #3rd level network
    if new_network.network.length > 1
      new_network.network.drop(1).each do |interactor_accession|
        new_network.add_interactions(interactor_accession) #Obtain interactors of the interactors
      end
    end
    
    new_network.protein_list = new_network.protein_list.sort_by(&:gene_id) #Sort the proteins from the list by the gene id
    new_network.add_filtered_network #Filter the network
  end  
end