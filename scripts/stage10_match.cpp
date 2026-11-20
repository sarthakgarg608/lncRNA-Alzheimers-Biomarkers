// Exhaustive ungapped 50-mer searches at Hamming distance <=2 on both strands.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
struct Seed { int query, offset; };
struct Hit { int distance; string transcript, strand, type; };
int base(char c) {
    switch(c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}
int main(int argc, char** argv) {
    if (argc != 2) return 1;
    ios::sync_with_stdio(false);
    ifstream file(argv[1]);
    vector<string> queries;
    unordered_multimap<uint32_t, Seed> seeds;
    // A small bit filter avoids most hash lookups when scanning the reference.
    vector<bool> present(1<<24);
    string seq;
    while (file >> seq) {
        if (seq.size()!=50 || seq.find_first_not_of("ACGT")!=string::npos) return 2;
        for (int orientation=0; orientation<2; ++orientation) {
            int q=queries.size();
            queries.push_back(seq);
            for (int offset : {0,17,34}) {
                uint32_t key=0;
                for (int j=0; j<16; ++j) key=(key<<2)|base(seq[offset+j]);
                seeds.emplace(key,Seed{q,offset});
                present[key & 0xffffff]=true;
            }
            reverse(seq.begin(),seq.end());
            for (char& c : seq) c="TGCA"[base(c)];
        }
    }
    if (queries.empty()) return 3;
    map<pair<int,string>,Hit> hits;
    string line, transcript, gene, type, reference;
    auto scan = [&]() {
        uint32_t key=0;
        int run=0;
        for (int i=0; i<(int)reference.size(); ++i) {
            int b=base(reference[i]);
            if (b<0) { run=0; key=0; continue; }
            key=(key<<2)|b;
            if (++run<16 || !present[key & 0xffffff]) continue;
            auto range=seeds.equal_range(key);
            for (auto it=range.first; it!=range.second; ++it) {
                auto seed=it->second;
                int start=i-15-seed.offset;
                if (start<0 || start+50>(int)reference.size()) continue;
                int distance=0;
                for (int j=0; j<50 && distance<=2; ++j)
                    distance+=(queries[seed.query][j]!=reference[start+j]);
                if (distance>2) continue;
                auto id=make_pair(seed.query/2,gene);
                auto old=hits.find(id);
                if (old==hits.end() || distance<old->second.distance)
                    hits[id]={distance,transcript,seed.query%2?"-":"+",type};
            }
        }
    };
    while (getline(cin,line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0]!='>') {reference+=line; continue;}
        scan(); reference.clear();
        stringstream header(line.substr(1));
        vector<string> fields; string field;
        while (getline(header,field,'|')) fields.push_back(field);
        if (fields.size()<8) return 4;
        transcript=fields[0]; gene=fields[1]; type=fields[7];
    }
    scan();
    for (auto const& row : hits) {
        auto const& h=row.second;
        cout<<queries[2*row.first.first]<<'\t'<<row.first.second<<'\t'<<h.distance<<'\t'
            <<h.transcript<<'\t'<<h.strand<<'\t'<<h.type<<'\n';
    }
    return cin.bad() || !cout ? 5 : 0;
}

