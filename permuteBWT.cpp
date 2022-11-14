#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <algorithm>
#include <cstring>
#include <cassert>

using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::ifstream;
using std::pair;
using std::make_pair;
using std::vector;

// ASCII alphabet size
const int alph_size = 128;

// empty all Stack and write the content in BWTblock vector
void permute_stack(std::stack< uint64_t* > &Stack, std::vector< pair <char,uint64_t> > &BWTblock){
    // extract the first element
    size_t i;
    int matchingInd;
    uint64_t* prev = Stack.top();
    Stack.pop();
    if( Stack.size() > 0 ){
        // iterate until the Stack is empty 
        while( !Stack.empty() ){
            // extract current Parick vector
            uint64_t* curr = Stack.top();
            Stack.pop();
            matchingInd = -1;
            // check for characters matching with previous 
            // Parick vector
            for(i=0; i<alph_size; ++i){
                if(prev[i] > 0){
                    if( curr[i] == 0 || (curr[i] > 0 && matchingInd !=-1 ) ){ 
                        // insert character if not matching
                        BWTblock.push_back(  make_pair(i, prev[i]) );
                        prev[i]=0;
                    }
                    else{ matchingInd=i; }
                }
            }
            // if there was a matching character
            if(matchingInd != -1){ 
                BWTblock.push_back( make_pair(matchingInd, prev[matchingInd]+curr[matchingInd]) ); 
                prev[matchingInd] = curr[matchingInd] = 0;
            } 
            // assing to the previous Parick vector
            // the current one
            delete[] prev;
            prev=curr; 
        } 
    }  
    // process the last Parick vector
    for(i=0; i<alph_size; ++i){ 
        if(prev[i] > 0){ BWTblock.push_back( make_pair(i,prev[i]) ); }
        prev[i]=0;
    }
    // free prev Parick vector
    delete[] prev;
}

// write the BWT runs in BWTblock vector
void write_block(std::vector< pair <char,uint64_t> > &BWTblock, int &last, FILE *obwt){
    // initialize variables
    int64_t i; uint64_t rep;
    char c;
    last = BWTblock[0].first;
    vector<uint8_t> run_buffer;
    // write all BWT runs
    for(i=BWTblock.size()-1; i>-1; --i){
    	c = BWTblock[i].first;
    	rep = BWTblock[i].second;
    	// write a run in one pass
    	if(run_buffer.size() < rep){ run_buffer.resize(rep); }
    	memset(&run_buffer[0],c,rep*sizeof(char));
	    if(fwrite(&run_buffer[0], sizeof(char), rep, obwt) != rep ){
	        std::cerr << "Error writing to optimal BWT file!\n"; exit(1);
	    }
	}
    // clean BWT segment vectors
    run_buffer.clear();
    BWTblock.clear();
}

void processP(uint64_t* Parikh, std::stack< uint64_t* > &Stack, std::vector< pair <char,uint64_t> > &BWTblock,
              int &last, FILE *obwt){
    // process the current Parikh vector that is
    // for sure an interesting interval
    if( Stack.empty() ){
        //std::cout << "Stack empty\n";
        // check last inserted char
        if( Parikh[last] > 0 ){
            // write a new BWT block
            BWTblock.push_back( make_pair( last, Parikh[last] ) );
            Parikh[ last ] = 0;
            write_block(BWTblock, last, obwt);
        }
        // insert a new Parikh vector in the stak
        uint64_t* newParikh = new uint64_t[alph_size]();
        memcpy(&newParikh[0],&Parikh[0],alph_size*sizeof(uint64_t));
        Stack.push(newParikh);
        last = -1;
    }
    else
    {
        // if the Stack is not empty
        uint64_t* top = Stack.top();
        int64_t m = 0, d = -1;
        // check number of matching characters
        for(int y = 0; y < alph_size; ++y){
            // if character j occurs at least once
            if( Parikh[y] > 0 ){ 
                // check if character is matching in the previous
                // Parick vector
                if( top[y] > 0 ){ m++; d = y;}
            }
        }
        // if we have less than two matching characters
        if( m < 2 ){
            // only one matching characters 
            if( m == 1 )
            {
                BWTblock.push_back( make_pair(d, top[d]+Parikh[d]) );
                top[d] = Parikh[d] = 0;
            }
            // permute the characters in the vector
            permute_stack(Stack,BWTblock);
            write_block(BWTblock, last, obwt);
        }
        // push a new Parikh vector
        uint64_t* newParikh = new uint64_t[alph_size]();
        memcpy( newParikh, Parikh, alph_size*sizeof(uint64_t) );
        Stack.push(newParikh);
    }
    // reset the Parikh vector
    memset(Parikh,0,alph_size*sizeof(uint64_t));
}

int main(int argc, char** argv)
{

	// show script usage information
	if(argc<4){
		cout << "\nUsage: " << argv[0] << " input.bwt input.sap buffer_size\n";
		cout << "Compute the optimal BWT of the input file and output it to input.optbwt\n";
		cout << "it takes in input the paths of the BWT and the SAP array and the size of the memory buffer (in MB)\n";
		exit(1);
	}
	else{
		puts("==== Command line:");
		for(int i=0;i<argc;i++)
			cout << argv[i] << " ";
		puts("");
	}

	// compute input buffer size
	size_t buffer_size = (atoll(argv[3])*1000000)/2;
	// initialize the input buffer
	uint8_t *bwt_buffer = new uint8_t[buffer_size];
	uint8_t *sap_buffer = new uint8_t[buffer_size];
	// read the input paths
	string bwt_file = argv[1];
	string sap_file = argv[2];
	// initialize input file streams
    ifstream fin_bwt(bwt_file);
    ifstream fin_sap(sap_file);

    // initialize variables
    int64_t i, zero_start = 0;
    int64_t count_bwt, count_sap;
    int last = 0, last_bwt = 0, last_sap = 48;
    // initialize stack data strcture
    std::stack< uint64_t* > Stack;
    // initialize current and previous Parick vector
    uint64_t* Parikh = new uint64_t[alph_size]();
    // initialize BWT block vector
    std::vector< pair <char,uint64_t> > BWTblock;
    // initalize output files
    std::string obwtfile = bwt_file + ".optbwt";
    FILE *obwt;
    // open output file
    if((obwt = fopen(obwtfile.c_str(), "w+")) == nullptr){
        std::cerr << "open() file " + obwtfile + " failed!" << std::endl;
        exit(1);
    }

    // start reading the input files from stream
    while(fin_bwt){
    	// read the next data chunk
		fin_bwt.read(reinterpret_cast<char*>(bwt_buffer), buffer_size);
		fin_sap.read(reinterpret_cast<char*>(sap_buffer), buffer_size);
	  	// Get the number of bytes actually read
	  	count_bwt = fin_bwt.gcount();
	  	count_sap = fin_sap.gcount();
	  	// if it has been read a different number of elements break
	  	if( count_bwt != count_sap ){ cerr << "Error! The input files have different lengths... exiting.\n"; exit(1); }
	  	// If nothing has been read, break
	  	if( !count_bwt ){ break; }
        // resent i index
        i = 0, zero_start = 0;
        // manage previous buffer border
        if( last_sap == 48 ){
            // check current SAP value
            if( sap_buffer[i] == 48 ){
                if( last_bwt > 0 ){
                    // check stack
                    if(!Stack.empty()){
                        uint64_t* top = Stack.top();
                        if( top[last_bwt] > 0 ){
                            BWTblock.push_back( make_pair( last_bwt, top[last_bwt] ) );
                            top[ last_bwt ] = 0;
                        }
                        permute_stack(Stack,BWTblock);
                        write_block(BWTblock, last, obwt);
                    }
                    putc(last_bwt, obwt); last = last_bwt; 
                }
            }
            else
            {  
                // add last bwt char to Parikh vector
                Parikh[last_bwt]++;
                // we are in a interesting interval
                while( i < count_bwt ){
                    // stop when seeing a zero
                    if(sap_buffer[i] == 48){ break; }
                    // increment Parikh counter
                    Parikh[bwt_buffer[i]]++;
                    // increment counter
                    i++;
                }
                if( i < count_bwt ){
                    // process last Parikh vector
                    processP(Parikh, Stack, BWTblock, last, obwt);
                    // update start
                    zero_start = i;
                }
            }
        }
        else{
            // last sap was 1
            if( sap_buffer[i] == 48 ){
                // process last Parikh vector
                processP(Parikh, Stack, BWTblock, last, obwt);
            }
            else
            {
                // we are in a interesting interval
                while( i < count_bwt ){
                    // stop when seeing a zero
                    if( sap_buffer[i] == 48 ){ break; }
                    // increment Parikh counter
                    Parikh[bwt_buffer[i]]++;
                    // increment counter
                    i++;
                }
                if( i < count_bwt ){
                    // process last Parikh vector
                    processP(Parikh, Stack, BWTblock, last, obwt);
                    // update start
                    zero_start = i;
                }
            }
        }
	  	// scan all positions of the buffer
        while( i < count_bwt ){
            // search next interesting interval
            while( i < count_bwt ){
                // stop when finding 1 value
                if( sap_buffer[i] == 49 ){ 
                    // increment Parikh vector counters
                    Parikh[bwt_buffer[i-1]]++;
                    break;
                }
                // increment counter
                i++;
            } 
            // check if there is something to write
            if( i-2 >= zero_start ){
                // if Stack is not empty
                if(!Stack.empty()){
                    uint64_t* top = Stack.top();
                    if(top[bwt_buffer[zero_start]] > 0){
                        BWTblock.push_back( make_pair( bwt_buffer[zero_start], top[bwt_buffer[zero_start]] ) );
                        top[ bwt_buffer[zero_start] ] = 0;
                    }
                    permute_stack(Stack,BWTblock);
                    write_block(BWTblock, last, obwt);
                }
                // write a sequence of zeros
                if( (int64_t)fwrite(&bwt_buffer[zero_start],sizeof(char),(i-1-zero_start),obwt) != (i-1-zero_start) ){
                    std::cerr << zero_start << " " << (i-1-zero_start) << " " << i << std::endl;
                    std::cerr << "Error writing in the optimalBWT file... exiting.\n";
                    exit(1);
                }
                // update last
                last = bwt_buffer[i - 2];
            } 
            // scan whole SAP-interval
            while( i < count_bwt ){
                // stop when finding 1 value
                if( sap_buffer[i] == 48 ){
                    // process the last Parikh vector
                    processP(Parikh, Stack, BWTblock, last, obwt);
                    zero_start = i;
                    break;
                }
                // increment Parikh vector counter
                Parikh[bwt_buffer[i]]++;
                // increment counter
                i++;
            }      
        }
        // set last bwt and sap values
        last_bwt = bwt_buffer[count_bwt-1];
        last_sap = sap_buffer[count_bwt-1];
    }
    // write the last Parick vector
    if( last_sap == 49 ){
        // process last Parikh vector
        processP(Parikh, Stack, BWTblock, last, obwt);
        //
        if(!Stack.empty()){
            permute_stack(Stack,BWTblock); 
            write_block(BWTblock, last, obwt);
        }
    }
    else
    {
        if(!Stack.empty()){
            uint64_t* top = Stack.top();
            if(top[ last_bwt ] > 0){
                BWTblock.push_back( make_pair( last_bwt, top[ last_bwt ] ) );
                top[ last_bwt ] = 0;
            }
            permute_stack(Stack,BWTblock);
            write_block(BWTblock, last, obwt);
        }
        // write last char
        putc(last_bwt, obwt);
    }
    // close output file
    fclose(obwt);
    delete[] bwt_buffer;
    delete[] sap_buffer;

    return 0;
}