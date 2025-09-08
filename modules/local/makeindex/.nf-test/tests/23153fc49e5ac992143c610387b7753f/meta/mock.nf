import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { MAKEINDEX } from '/data/ABO-analysis-31-Oct-2023/nf-core-abotyper/modules/local/makeindex/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                input[0] = [
                    [ id:'test.complex.name.fasta', exon:'exon6' ],
                    file(params.pipelines_testdata_base_path + '/abotyper/abotyper/refs/A1_01_01_1_reference_Exon6.fasta.fai', checkIfExists: true)
                ]
                input[1] = [
                    [ id:'another.test.fasta', exon:'exon7' ],
                    file(params.pipelines_testdata_base_path + '/abotyper/abotyper/refs/A1_01_01_1_reference_Exon7.fasta.fai', checkIfExists: true)
                ]
                
    //----

    //run process
    MAKEINDEX(*input)

    if (MAKEINDEX.output){

        // consumes all named output channels and stores items in a json file
        for (def name in MAKEINDEX.out.getNames()) {
            serializeChannel(name, MAKEINDEX.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = MAKEINDEX.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
