#include "read_geom.h"


void read_geom_line(vector<Atom> &Atoms, string line){
    Atom temp;
    if (line.find("Nuclear Charges and Cartesian Coordinates:", 0)!= string::npos) return;
    if (line.find("---------------------------", 0)!= string::npos) return;
    if (line.find("Atom      Charge                X                  Y                  Z", 0)!= string::npos) return;
    if (line == " ") return;
    std::istringstream iss(line);
    iss >> temp.ELEMENT >>temp.charge>> temp.pos.x >> temp.pos.y >> temp.pos.z;;
    Atoms.push_back(temp);
}


void read_geometry(vector<Atom> &Atoms, string filename){
    ifstream inp;
    string line;
    bool top = false;
    bool buttom = false;
    inp.open(filename.c_str());
    while(getline(inp, line)){
        if (line.find("Nuclear Charges and Cartesian Coordinates", 0) != string::npos) top = true;
        if (line.find("Some Atomic Properties:", 0) != string::npos) top = false;
        if (top){
            read_geom_line(Atoms, line);
        }
    }

}


void read_FluxJ_line(vector<FluxJ> &FluxJs, string line){
    FluxJ temp;
    if (line.find("                     xxxx                    ", 0)!= string::npos) return;

    std::istringstream iss(line);
    iss >> temp.IELEMENT >>temp.IIELEMENT>> temp.J.x >> temp.J.y >> temp.J.z;;
    FluxJs.push_back(temp);
}

void read_extracted_FluxJ(vector<FluxJ> &FluxJs, string filename){
    ifstream inp;
    string line;
    inp.open(filename.c_str());
    while(getline(inp, line)) read_FluxJ_line(FluxJs, line);
}