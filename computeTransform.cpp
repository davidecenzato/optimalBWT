#include "computeTransform.h"

// empty all Stack and write the content in BWTblock vector
void permute_stack(std::stack< uint_t* > &Stack, std::vector< std::pair <char,uint_t> > &BWTblock){
    // extract the first element
    size_t i;
    int matchingInd;
    uint_t* prev = Stack.top();
    Stack.pop();
    if( Stack.size() > 0 ){
        // iterate until the Stack is empty 
        while( !Stack.empty() ){
            // extract current Parick vector
            uint_t* curr = Stack.top();
            Stack.pop();
            matchingInd = -1;
            // check for characters matching with previous 
            // Parick vector
            for(i=0; i<alph_size; ++i){
                if(prev[i] > 0){
                    if( curr[i] == 0 || (curr[i] > 0 && matchingInd !=-1 ) ){ 
                        // insert character if not matching
                        BWTblock.push_back(  std::make_pair(i, prev[i]) );
                        prev[i]=0;
                    }
                    else{ matchingInd=i; }
                }
            }
            // if there was a matching character
            if(matchingInd != -1){ 
                BWTblock.push_back( std::make_pair(matchingInd, prev[matchingInd]+curr[matchingInd]) ); 
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
        if(prev[i] > 0){ BWTblock.push_back( std::make_pair(i,prev[i]) ); }
        prev[i]=0;
    }
    // free prev Parick vector
    delete[] prev;
}

// write the BWT runs in BWTblock vector
void write_block(std::vector< std::pair <char,uint_t> > &BWTblock, uint8_t* BWT, uint_t &cnt, int &last){
    // initialize variables
    int64_t i, rep;
    char c;
    last = BWTblock[0].first;
    // write all BWT runs
    for(i=BWTblock.size()-1; i>-1; --i){
        c = BWTblock[i].first;
        rep = BWTblock[i].second;
        // write a BWT run
        memset(&BWT[cnt], c, rep*sizeof(char));
        cnt += rep;
    }
    // clean BWT segment vector
    BWTblock.clear();
}

void permute_bwt_unsigned(std::vector<char> &Text, std::vector<uint_t> &SA, std::vector<uint8_t> &I,
                          std::string basename, uint_t n){

    // initialize stack data strcture
    std::stack< uint_t* > Stack;
    // initalize output files
    std::string bwtfile = basename + std::string(".optbwt");
    FILE *obwt;
    // initialize variables
    size_t i, y;
    int j = 0, m, d;
    int last = 0, no_char = 0; 
    // initialize BWT pointer to SA
    uint_t cnt = 0;
    uint8_t* ptr = (uint8_t*)&SA[cnt];
    // initialize current and previous Parick vector
    uint_t* Parikh = new uint_t[alph_size]();
    // initialize BWT segment vector
    std::vector< std::pair <char,uint_t> > BWTblock;
    // insert eh first interval
    // store current BWT char
    j = (SA[0] > 0)?Text[SA[0]-1]:'$';  
    if( j==1 ){ j='$'; }
    no_char += (Parikh[j]++)==0?1:0;
    // increment the current SAP interval 
    // scan the suffix array while identifying the 
    // interesting intervals
    for(i = 1; i < n; ++i){
        // if we are in at the beginning SAP interval
        if( iget(i) == 0 ){
            // if the interval is not interesting
            if( no_char == 1 ){
                // update last character
                last = j;
                // if the stack is empty
                if( Stack.empty() ){
                    memset(&ptr[cnt], j, Parikh[j]*sizeof(char));
                    cnt += Parikh[j]; Parikh[j] = 0;
                }
                else
                {
                    uint_t* top = Stack.top();
                    // if there is a matching char with previous interval
                    BWTblock.push_back( std::make_pair(j,top[j]+Parikh[j]) );
                    top[j]=0; Parikh[j]=0;
                    // empty the stack
                    permute_stack( Stack, BWTblock );
                    write_block( BWTblock, ptr, cnt, last );
                }
                // insert current char and continue
                j = (SA[i] > 0)?Text[SA[i]-1]:'$'; 
                if( j==1 ){ j='$'; }
                Parikh[j] = 1; no_char = 1;
                continue;
            }
            else if( no_char > 1 )
            {
                // we are in an interesting interval
                // when the Stack is empty
                if( Stack.empty() )
                { 
                    assert(last != -1);
                    // if there is a matching character
                    if( Parikh[last] > 0 ){
                        // write a new BWT block
                        memset(&ptr[cnt], last, Parikh[last]*sizeof(char));
                        cnt += Parikh[last];
                        Parikh[ last ] = 0;
                    }
                    // insert a new Parikh vector in the stak
                    uint_t* newParikh = new uint_t[alph_size]();
                    memcpy(&newParikh[0],&Parikh[0],alph_size*sizeof(uint_t));
                    Stack.push(newParikh);
                    last = -1;
                }
                else
                {
                    // if the Stack is not empty
                    uint_t* top = Stack.top();
                    m = 0, d = -1;
                    // check number of matching characters
                    for(y = 0; y < alph_size; ++y){
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
                            BWTblock.push_back( std::make_pair(d, top[d]+Parikh[d]) );
                            top[d] = Parikh[d] = 0;
                        }
                        // permute the characters in the vector
                        permute_stack(Stack,BWTblock);
                        write_block( BWTblock, ptr, cnt, last );
                    }
                    // push a new Parikh vector
                    uint_t* newParikh = new uint_t[alph_size]();
                    memcpy( newParikh, Parikh, alph_size*sizeof(uint_t) );
                    Stack.push(newParikh);
                }

            }
            // reset the Parikh vector
            memset(Parikh,0,alph_size*sizeof(uint_t));
            no_char = 0;
        }
        // store current BWT char
        j = (SA[i] > 0)?Text[SA[i]-1]:'$';  
        if( j==1 ){ j='$'; }
        // increment the current SAP interval
        no_char += (Parikh[j]++)==0?1:0; 
    }
    // write the last Parick vector
    Stack.push(Parikh);
    permute_stack(Stack, BWTblock);
    write_block( BWTblock, ptr, cnt, last );
    // open output file
    if((obwt = fopen(bwtfile.c_str(), "w+")) == nullptr){
        std::cerr << "open() file " + bwtfile + " failed!" << std::endl;
        exit(-1);
    }
    // write the optimal BWT
    if(fwrite(ptr, sizeof(char), n, obwt)!=n){ 
        std::cerr << "Error writing in the optimal BWT file!" << std::endl; 
        exit(1); }
        
    // close optimal bwt file
    fclose(obwt);
}

// write the BWT and SAP array to file
void write_bwt_sap(std::vector<char> &Text, std::vector<uint_t> &SA, std::vector<uint8_t> &I,
                   std::string basename, uint_t n){

    // intialize variables
    char j, b;
    // initalize output files
    std::string bwtfile, sapfile;
    FILE *bwt, *sap;

    uint8_t* ptr = (uint8_t*)&SA[0];

    bwtfile = basename + std::string(".mdollarbwt");
    sapfile = basename + std::string(".sap");
    // open output files
    if((bwt = fopen(bwtfile.c_str(), "w+")) == nullptr){
        std::cerr << "open() file " + bwtfile + " failed!" << std::endl;
        exit(-1);
    }
    if((sap = fopen(sapfile.c_str(), "w+")) == nullptr){
        std::cerr << "open() file " + bwtfile + " failed!" << std::endl;
        exit(-1);
    }
    // scan the suffix array and write the BWT
    for(uint_t i = 0; i < n; ++i) {
        // compute the new BWT char
        j = (SA[i] > 0)?Text[SA[i]-1]:'$';
        if( j==1 ){ j='$'; }
        ptr[i] = j;
        b = iget(i);
        putc(b, sap); 
    }
    if(fwrite(ptr, sizeof(char), n, bwt)!=n){ std::cerr << "grossi errori" << std::endl; exit(1); }
    // close output files
    fclose(bwt);
    fclose(sap);
}

// write the BWT to file
void write_bwt(std::vector<char> &Text, std::vector<uint_t> &SA, std::string basename, uint_t n){

    // intialize variables
    char j, b;
    // initalize output files
    std::string bwtfile, sapfile;
    FILE *bwt;

    uint8_t* ptr = (uint8_t*)&SA[0];

    bwtfile = basename + std::string(".bwt");
    // open output files
    if((bwt = fopen(bwtfile.c_str(), "w+")) == nullptr){
        std::cerr << "open() file " + bwtfile + " failed!" << std::endl;
        exit(-1);
    }
    // scan the suffix array and write the BWT
    for(uint_t i = 0; i < n; ++i) {
        // compute the new BWT char
        j = (SA[i] > 0)?Text[SA[i]-1]:'$';
        if( j==1 ){ j='$'; }
        ptr[i] = j;
    }
    if(fwrite(ptr, sizeof(char), n, bwt)!=n){ std::cerr << "Error writing in the BWT file!" << std::endl; exit(1); }
    // close output files
    fclose(bwt);
}

// function for constructing and storing the optimalBWT of a string collection
// light version
void compute_optbwt_unsigned(Args arg){ 
    // initialize starting time
    auto start = std::chrono::steady_clock::now();
    // input text initialization
    std::vector<char> Text;
    //std::stack< std::array<uint_s,alph_size>* > Stack;
    uint_t n=0, ns=0;
    // initialize output basename
    std::string basename = arg.outname;
    // initialize boolean values
    int format = arg.format;
    ////bool printsa = arg.sa;
    // read input
    switch(format) {
        case 0:
            //std::cout << "load fasta file" << std::endl;
            // load a fasta file
            load_fasta_conc(arg.filename.c_str(),Text,n,ns);
            break;
        case 1:
            // load a fastq file
            load_fastq_conc(arg.filename.c_str(),Text,n,ns);
            break;
        default:
            std::cerr << "error! please select a valid format." << std::endl;
            exit(-1);
    }
    // initialize the circular Suffix array
    std::vector<uint_t> SA(n,0);
    // initialize the SAP-array
    std::vector<uint8_t> I((n/8)+1,0);
    // compute the suffix array
    optsais_unsigned_sa_sap(&Text[0], &SA[0], &I[0], n, ns, alph_size, 1);

    auto end = std::chrono::steady_clock::now();
    std::cout << "Constructing the SA: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;
    start = std::chrono::steady_clock::now();

    // permute the BWT characters using the Bentley et al. algorithm
    permute_bwt_unsigned(Text, SA, I, basename, n);

    end = std::chrono::steady_clock::now();
    std::cout << "Permuting the characters: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;

    // free memory
    Text.clear();
    SA.clear();
    I.clear();
}


// function for constructing and storing the optimalBWT of a string collection
// light version
void compute_bwt_sap_unsigned(Args arg){ 
    // initialize starting time
    auto start = std::chrono::steady_clock::now();
    // input text initialization
    std::vector<char> Text;
    //std::stack< std::array<uint_s,alph_size>* > Stack;
    uint_t n=0, ns=0;
    // initialize output basename
    std::string basename = arg.outname;
    // initialize boolean values
    int format = arg.format;
    ////bool printsa = arg.sa;
    // read input
    switch(format) {
        case 0:
            //std::cout << "load fasta file" << std::endl;
            // load a fasta file
            load_fasta_conc(arg.filename.c_str(),Text,n,ns);
            break;
        case 1:
            // load a fastq file
            load_fastq_conc(arg.filename.c_str(),Text,n,ns);
            break;
        default:
            std::cerr << "error! please select a valid format." << std::endl;
            exit(-1);
    }
    // initialize the circular Suffix array
    std::vector<uint_t> SA(n,0);
    // initialize the SAP-array
    std::vector<uint8_t> I((n/8)+1,0);
    // compute the suffix array
    optsais_unsigned_sa_sap(&Text[0], &SA[0], &I[0], n, ns, alph_size, 1);
    
    auto end = std::chrono::steady_clock::now();
    std::cout << "Constructing the SA: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;

    // permute the BWT characters using the Bentley et al. algorithm
    start = std::chrono::steady_clock::now();
    write_bwt_sap(Text, SA, I, basename, n);
    end = std::chrono::steady_clock::now();
    std::cout << "Writing BWT and SAP: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;

    // free memory
    Text.clear();
    SA.clear();
    I.clear();
}

// function for constructing and storing the optimalBWT of a string collection
// light version
void compute_bwt_unsigned(Args arg){ 
    // initialize starting time
    auto start = std::chrono::steady_clock::now();
    // input text initialization
    std::vector<char> Text;
    //std::stack< std::array<uint_s,alph_size>* > Stack;
    uint_t n=0, ns=0;
    // initialize output basename
    std::string basename = arg.outname;
    // initialize boolean values
    int format = arg.format;
    ////bool printsa = arg.sa;
    // read input
    switch(format) {
        case 0:
            // load a fasta file
            load_fasta_conc(arg.filename.c_str(),Text,n,ns);
            break;
        case 1:
            // load a fastq file
            load_fastq_conc(arg.filename.c_str(),Text,n,ns);
            break;
        default:
            std::cerr << "error! please select a valid format." << std::endl;
            exit(-1);
    }
    // initialize the circular Suffix array
    std::vector<uint_t> SA(n,0);
    // compute the suffix array
    optsais_unsigned_sa(&Text[0], &SA[0], n, ns, alph_size, 1);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Constructing the SA: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;

    // permute the BWT characters using the Bentley et al. algorithm
    start = std::chrono::steady_clock::now();
    write_bwt(Text, SA, basename, n);
    end = std::chrono::steady_clock::now();
    std::cout << "Writing BWT: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;

    // free memory
    Text.clear();
    SA.clear();
}