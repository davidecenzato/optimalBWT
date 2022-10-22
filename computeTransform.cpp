#include "computeTransform.h"

// empty all Stack and write the content in BWTseg vector
void empty_stack_unsigned(std::stack< uint_t* > &Stack, std::vector<char> &BWTseg){
    // extract the first element
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
            for(size_t x=0; x<alph_size; ++x){
                if(prev[x] > 0){
                    if( curr[x] == 0 || (curr[x] > 0 && matchingInd !=-1 ) ){ 
                        // insert character if not matching
                        BWTseg.insert(BWTseg.end(), prev[x], x);
                        prev[x]=0;
                    }
                    else{ matchingInd=x; }
                }
            }
            // if there was a matching character
            if(matchingInd != -1){ 
                BWTseg.insert(BWTseg.end(), prev[matchingInd]+curr[matchingInd], matchingInd); 
                prev[matchingInd] = curr[matchingInd] = 0;
            }
            // assing to the previous Parick vector
            // the current one
            //free(prev); 
            delete[] prev;
            prev=curr; 
        } 
    }  
    // process the last Parick vector
    for(uint_t x=0; x<alph_size; ++x){ 
        if(prev[x] > 0) BWTseg.insert(BWTseg.end(), prev[x], x);
        prev[x]=0;
    }
    // free prev Parick vector
    //free(prev);
    delete[] prev;
}

// write BWTseg
void write_segment(std::vector<char> &BWTseg, bool rev, int &last, FILE *obwt){
    // write BWT segment
    if(rev) std::reverse(BWTseg.begin(), BWTseg.end());
    last = BWTseg[BWTseg.size() - 1];
    if(fwrite(&BWTseg[0], sizeof(char), BWTseg.size(), obwt) != BWTseg.size()){
        std::cerr << "Error writing to optimal BWT file!\n"; exit(-1);
    }
    // clean BWT segment vector
    BWTseg.clear();
}

void permute_bwt_usigned(std::vector<char> &Text, std::vector<uint_t> &SA, std::vector<uint8_t> &I,
                 std::string basename, uint_t n){

    // initialize stack data strcture
    std::stack< uint_t* > Stack;
    // initalize output files
    std::string bwtfile;
    FILE *obwt;

    bwtfile = basename + std::string(".optbwt");
    // open output files
    if((obwt = fopen(bwtfile.c_str(), "w+")) == nullptr){
        std::cerr << "open() file " + bwtfile + " failed!" << std::endl;
        exit(-1);
    }
    // initialize current and previous Parick vector
    uint_t* Parick; 
    // initialize BWT segment vector
    std::vector<char> BWTseg;
    // char len 1 intervals
    int len1 = -1, last1 = alph_size-1;
    // scan the suffix array while identifying the 
    // interesting intervals
    int nozero = 0;
    for(uint_t i = 0; i < n-1; ++i) {
        // compute the new BWT char
        int j = (SA[i] > 0)?Text[SA[i]-1]:'$';  
        if( j==1 ){ j='$'; }
        // check if it is a length one interval or the starting
        // position of a new interval
        if(!iget(i)){
            // len one interval
            if(!iget(i+1)){
                // char contained in the len one interval
                len1 = j; last1 = j;
                // if stack is empty write the character
                if(Stack.empty()){ putc(j, obwt);  continue; }
                else{
                    uint_t* top = Stack.top();
                    // if there is a matching char with previous interval
                    BWTseg.insert(BWTseg.end(), top[len1]+1, len1);
                    top[len1]=0;
                    // empty the stack
                    empty_stack_unsigned(Stack,BWTseg);
                    write_segment(BWTseg,true,last1,obwt);
                }
            }
            // starting point of a new interval
            else{
                // initialize a new Parick vector
                Parick = new uint_t[alph_size]();
                // add char to Parick vector
                Parick[j] += 1;
                nozero = 1;
                // continue to the next iteration
                continue;
            }
        // I[i] = 1
        }else{
            // add char to Parick vector
            if(Parick[j]==0){ nozero++; }
            Parick[j] += 1;
            // interval elongation
            if(iget(i+1)){
                continue; 
            }
            // last position of the current interval
            else{
                assert(nozero > 0); 
                if(!Stack.empty()){
                    if(nozero > 1){
                        // take Parick vector at the top of stack
                        uint_t* top = Stack.top();
                        int matched = 0, matchingInd = -1;
                        // check number of matching characters
                        for(size_t y=0; y<alph_size; ++y){
                            // if character j occurs at least once
                            if(Parick[y]>0){ 
                                // check if character is matching in the previous
                                // Parick vector
                                if( top[y]>0 ){ matched++; matchingInd = y;}
                            }
                        }
                        // if there is exactly one match
                        if(matched == 1){
                            // add the the BWT the part matching with previous interval
                            BWTseg.insert(BWTseg.end(), top[matchingInd]+Parick[matchingInd], matchingInd);
                            top[matchingInd] = Parick[matchingInd] = 0;
                            // empty the stack
                            empty_stack_unsigned(Stack,BWTseg);
                            write_segment(BWTseg,true,last1,obwt);
                        }
                        else if (matched == 0){
                            // no matching characters empty the stack
                            empty_stack_unsigned(Stack,BWTseg);
                            if(BWTseg.size() > 0) write_segment(BWTseg,true,last1,obwt);
                        } 
                        // push the last Parick vector
                        Stack.push(Parick);
                        // push the current parick vector to the stack
                    }else{
                        // not interesting interval (only one character)
                        uint_t* top = Stack.top();
                        // if the character is matching with the previous interval
                        if(top[j] > 0){ 
                            // insert matching part
                            BWTseg.insert(BWTseg.end(), top[j], j);
                            top[j] = 0;
                            // empty stack
                            empty_stack_unsigned(Stack,BWTseg);
                            write_segment(BWTseg,true,last1,obwt);
                        }
                        Stack.push(Parick);
                    }
                // if the Stack is already empty
                }else{
                    assert(last1 != -1);
                    // if there is a matching character and the stack is empty
                    if(Parick[last1] > 0 && nozero > 1){ 
                        BWTseg.insert(BWTseg.end(), Parick[last1], last1);
                        Parick[last1] = 0;
                        write_segment(BWTseg, false, last1, obwt);
                    }
                    Stack.push(Parick); last1 = -1; 
                }
            }
        }
    }
    // process last suffix
    int j = (SA[n-1] > 0)?Text[SA[n-1]-1]:'$';  /// 011110011001 011110110010
    if( j==1 ){ j='$'; }
    // add last BWT segment
    if(iget(n-1)){
        // add char to Parick vector
        Parick[j] += 1;
        Stack.push(Parick);
        empty_stack_unsigned(Stack,BWTseg);
        write_segment(BWTseg,true,last1,obwt);
    }else{
        len1 = j;
        if(Stack.empty()){ putc(j, obwt); }
        else{
            uint_t* top = Stack.top();
            // if there is a matching char
            BWTseg.insert(BWTseg.end(), top[len1]+1, len1);
            top[len1]=0;
            empty_stack_unsigned(Stack,BWTseg);
            write_segment(BWTseg,true,last1,obwt);
        }
    }
    // close optimal bwt file
    fclose(obwt);

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
    optsais_unsigned(&Text[0], &SA[0], &I[0], n, ns, alph_size, 1);

    auto end = std::chrono::steady_clock::now();
    std::cout << "Constructing the SA: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;
    start = std::chrono::steady_clock::now();

    // permute the BWT characters using the Bentley et al. algorithm
    permute_bwt_usigned(Text, SA, I, basename, n);

    end = std::chrono::steady_clock::now();
    std::cout << "Permuting the characters: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " milliseconds" << std::endl;

    // free memory
    Text.clear();
    SA.clear();
    I.clear();
}

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
    optsais_unsigned(&Text[0], &SA[0], &I[0], n, ns, alph_size, 1);
    
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