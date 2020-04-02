#include "bio.h"
#include<iostream>
#include <string>
#include <vector>
using std::string;
    //project 2 by Morgan Mundell
    //CSE 232 section 010
    //genetics project

bool IsValidDNASequence(const std::string & input){
    int len = input.length();
    bool validity = true;

    for (int i = 0; i < len; i++){
        if (input.at(i) == 'A' || input.at(i) == 'C' || input.at(i) == 'T' || input.at(i) == 'G'){
            continue;
        }
        else {
            validity = false;
            break;
                //breaks the for loop if one of the chars in the string is not A C T or G
        }
    }
    return validity;
}

void GetReverseComplementSequence(const std::string & input,  std::string * const output){
    std::string seq_holder;
        //a string that can be manipulated (hence being called holder)
    const int len = static_cast<int> (input.length());

    for (int i = len; i > 0; i--){
        seq_holder.push_back(input.at(i - 1));

        if(seq_holder.at(len-i) == 'A'){
            seq_holder.at(len-i) = 'T';
            continue;
        }

        else if(seq_holder.at(len-i) == 'T'){
            seq_holder.at(len-i) = 'A';
            continue;
        }

        else if(seq_holder.at(len-i) == 'C'){
            seq_holder.at(len-i) = 'G';
            continue;
        }

        else if(seq_holder.at(len-i) == 'G'){
            seq_holder.at(len-i) = 'C';
            continue;
        }
    }

    *output = seq_holder;

}

void TUSwitch (std::string & input){
    //this function switches T and U, it is used for getting RNA transcript and getting the antiparallel strand

    const int len = static_cast<int> (input.length());

    for(int i = 0; i < len; i++){
            if (input.at(i) == 'T'){
                input.at(i) = 'U';
            }
    }
}

std::string GetRNATranscript(const std::string & input){
    std::string output;
    std::string * p_output = &output;

    GetReverseComplementSequence(input, p_output);

    TUSwitch(output);
        //after getting the reverse compliment the T's and U's must be switched
    return output;
}

std::vector<std::vector<std::string> > GetReadingFramesAsCodons(const std::string & input){
    std::vector<std::vector<std::string> > output;
    std::vector<std::string> v_holder;
        //hold strings in a vector to later be pushed back into the output vector
    const int len = static_cast<int> (input.length());
    std::string org_rna, antipara = input;
        //strings for the original RNA sequence and the antiparallel sequence

    org_rna = GetRNATranscript(input);

    TUSwitch(antipara);
        //antipara is the input with T and U switched

    for (int i = 0; i < 3; i++){
            v_holder.clear();
            int l = 0;

        for (int j = len - i; j >= 3; j -= 3){
            v_holder.push_back(org_rna.substr(i + 3*l, 3));

            l++;
                //l is used to count the times the inner for loop is run
        }
            output.push_back(v_holder);
    }           //this code is for the original RNA string

    for (int i = 0; i < 3; i++){
            v_holder.clear();
            int l = 0;
        for (int j = len - i; j >= 3; j -= 3){
            v_holder.push_back(antipara.substr(i + 3*l, 3));

            l++;
        }
            output.push_back(v_holder);
    }
            //this code is for the antiparallel string
    return output;
}

std::string Translate(const std::vector<std::string> & codon_sequence){
    std::string output;
    std::vector<std::string> c = codon_sequence;
    const int len = static_cast<int> (codon_sequence.size());

    for (int i = 0; i < len; i ++){
        if (c.at(i) == "GCU" || c.at(i) == "GCC" || c.at(i) == "GCA" || c.at(i) == "GCG"){
            output.push_back('A');
        }
        else if (c.at(i) == "CGU" || c.at(i) == "CGC" || c.at(i) == "CGA" || c.at(i) == "CGG" || c.at(i) == "AGA" || c.at(i) == "AGG"){
            output.push_back('R');
        }
        else if (c.at(i) == "AAU" || c.at(i) == "AAC"){
            output.push_back('N');
        }
        else if (c.at(i) == "GAU" || c.at(i) == "GAC"){
            output.push_back('D');
        }
        else if (c.at(i) == "UGU" || c.at(i) == "UGC"){
            output.push_back('C');
        }
        else if (c.at(i) == "CAA" || c.at(i) == "CAG"){
            output.push_back('Q');
        }
        else if (c.at(i) == "GAA" || c.at(i) ==  "GAG"){
            output.push_back('E');
        }
        else if (c.at(i) == "GGA" || c.at(i) ==  "GGC" || c.at(i) ==  "GGU" || c.at(i) ==  "GGG"){
            output.push_back('G');
        }
        else if (c.at(i) == "CAU" || c.at(i) ==  "CAC"){
            output.push_back('H');
        }
        else if (c.at(i) == "AUC" || c.at(i) ==  "AUA" || c.at(i) ==  "AUU"){
            output.push_back('I');
        }
        else if (c.at(i) == "CUU" || c.at(i) ==  "CUC" || c.at(i) ==  "UUA" || c.at(i) ==  "UUG" || c.at(i) ==  "CUA" || c.at(i) ==  "CUG"){
            output.push_back('L');
        }
        else if (c.at(i) == "AAA" || c.at(i) ==  "AAG"){
            output.push_back('K');
        }
        else if (c.at(i) == "AUG"){
            output.push_back('M');
        }
        else if (c.at(i) == "UUU" || c.at(i) ==  "UUC"){
            output.push_back('F');
        }
        else if (c.at(i) == "CCG" || c.at(i) ==  "CCA" || c.at(i) ==  "CCU" || c.at(i) ==  "CCC"){
            output.push_back('P');
        }
        else if (c.at(i) == "AGU" || c.at(i) ==  "AGC" || c.at(i) ==  "UCU" || c.at(i) ==  "UCC" || c.at(i) ==  "UCA" || c.at(i) ==  "UCG"){
            output.push_back('S');
        }
        else if (c.at(i) == "ACG" || c.at(i) ==  "ACA" || c.at(i) ==  "ACU" || c.at(i) ==  "ACC"){
            output.push_back('T');
        }
        else if (c.at(i) == "UGG"){
            output.push_back('W');
        }
        else if (c.at(i) == "UAU" || c.at(i) == "UAC"){
            output.push_back('Y');
        }
        else if (c.at(i) == "GUA" || c.at(i) == "GUG" || c.at(i) ==  "GUU" || c.at(i) ==  "GUC"){
            output.push_back('V');
        }
        else if (c.at(i) ==  "UAG" || c.at(i) ==  "UGA" || c.at(i) ==  "UAA"){
            output.push_back('*');
        }
            //info found in bio.h
            //statements to determine the amino acid strand
    }

    return output;
}

std::string GetLongestOpenReadingFrame(const std::string & DNA_sequence){
    int len = 0, pos = 0, max_len = 0;
    std::string dna;
    std::string long_str;
    std::vector<std::vector<std::string> > codons;

    codons = GetReadingFramesAsCodons(DNA_sequence);
        //converts DNA seq into condons, to later be used to find longest open reading frame
    int c_size = codons.size();

        for (int i = 0; i < c_size; i++){
            dna = Translate(codons.at(i));
            pos = dna.find("M");
                // info for .find() function found at:
                // http://www.cplusplus.com/reference/string/string/find/

            dna.erase(dna.begin(), dna.begin() + pos - 1);
            len = dna.find('*');
                // len is length of the reading frame
                // erased strand before opening 'M' to easily find the beginning

            if (len > max_len&& dna.find(" ") == std::string::npos){
                max_len = len;
                long_str = dna.substr(1, len);
                    //holds the longest open reading frame
            }
            dna.erase(dna.begin(), dna.begin() + len - 1);
                //erase so program does not find the same opening 'M'
            dna.clear();
        }

        std::cout << long_str;

    return long_str;
}

