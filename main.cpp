#include "read_geom.h"

struct d_flux_J{
    string IELEMENT;
    string IIELEMENT;
    vec prepJ; 
    double Jprep;
    
    vector<double> J2prep;


    void print(ostream& stream) const{
		stream<<this->IELEMENT<<"\t"<<this->IIELEMENT<<"\t";
		stream<<setprecision(11) << fixed << right << setw(20) <<this->Jprep << endl;
	}

    void print2(ostream& stream) const{
		stream<<this->IELEMENT<<"\t"<<this->IIELEMENT<<"\t";
        for(int i=0; i< this->J2prep.size(); ++i){
            stream<<setprecision(1) << fixed << right << setw(20) <<this->J2prep[i] << "\t";    
        }
		stream << endl;
	}

};


vec pos_element(string &element, vector<Atom> &Atoms){
    vec temp;
    for(int i=0; i<Atoms.size();++i){
        if (element == Atoms[i].ELEMENT) return Atoms[i].pos;
    }
    cout<<element<< "Warning - no atom found"<<endl;
    return temp;
}

vector<string> connected_atoms(string &element, vector<FluxJ> &FluxJs){
    vector<string> res;
    for (int i =0; i<FluxJs.size(); ++i){
        if(FluxJs[i].IELEMENT == element) res.push_back(FluxJs[i].IIELEMENT);
    }
    return res;
}

vec prep_unit_vector_on_atom(string ELEMENT, vector<Atom> &Atoms, vector<string> &con_atoms){
    vector<vec> vec_con_atoms_bonds;
    for(int i=0; i<con_atoms.size(); ++i){
        vec_con_atoms_bonds.push_back(diff_vec(pos_element(con_atoms[i], Atoms), pos_element(ELEMENT, Atoms)));
    }
    
    vector<vec> vec_prep_bonds;

    for(int i=0; i<vec_con_atoms_bonds.size(); ++i){
        for(int j =i+1; j<vec_con_atoms_bonds.size(); ++j){
            vec tempvec;
            tempvec = unit_vec(crossproduct(
                vec_con_atoms_bonds[i],
                vec_con_atoms_bonds[j]));
            if(tempvec.z >= 0) {
                vec_prep_bonds.push_back(tempvec);
                } else {
                    vec_prep_bonds.push_back(scaleproduct(-1.0, tempvec));
                }
        }
    }
    return unit_vec(average_vectors(vec_prep_bonds));
}


vector<vec> prep_unit_vector(vector<vector<string>> &ref_atoms, vector<Atom> &Atoms, vector<FluxJ> FluxJs){
    vector <vec> des_unit_vectors;
    for (int i=0; i<  ref_atoms.size(); ++i){
        vector<vec> prep_vecs;
        for(int j=0; j<ref_atoms[i].size(); ++j){
            vector<string> con_atoms;
            con_atoms = connected_atoms(ref_atoms[i][j], FluxJs);
            prep_vecs.push_back(prep_unit_vector_on_atom(ref_atoms[i][j], Atoms, con_atoms));
            con_atoms.clear();
        }
        des_unit_vectors.push_back(unit_vec(average_vectors(prep_vecs)));
        prep_vecs.clear();
    }
    return des_unit_vectors;
}




void calculate_desired_vectors(d_flux_J &temp, vector<Atom> &Atoms, FluxJ &FJ, vector<string> &Icon_atoms, vector<string> &IIcon_atoms){
//    vector<vec> vec_con_atoms_bonds;
//    vec bond = diff_vec(pos_element(FJ.IIELEMENT, Atoms), pos_element(FJ.IELEMENT, Atoms));
//    temp.bdirJ = projvecJonR(FJ.J,bond);
//    temp.Jbdir=veclength(temp.bdirJ);
//    if(veclength(sum_vec(bond, temp.bdirJ)) < veclength(bond)) temp.Jbdir = temp.Jbdir*-1.0;
/*    for(int i=0; i<Icon_atoms.size(); ++i){
        vec_con_atoms_bonds.push_back(diff_vec(pos_element(Icon_atoms[i], Atoms), pos_element(FJ.IELEMENT, Atoms)));
    }

    vector<vec> J_vec_prep_bonds;
    
    for(int i=0; i<vec_con_atoms_bonds.size(); ++i){
        for(int j =i+1; j<vec_con_atoms_bonds.size(); ++j){
            vec tempvec;
            tempvec = unit_vec(projvecJonR(FJ.J, crossproduct(
                vec_con_atoms_bonds[i],
                vec_con_atoms_bonds[j])));
            if(tempvec.z >= 0) {
                J_vec_prep_bonds.push_back(tempvec);
                } else {
                    J_vec_prep_bonds.push_back(scaleproduct(-1.0, tempvec));
                }
        }
    }
    temp.prepJ=unit_vec(average_vectors(J_vec_prep_bonds));
 */   //if (temp.prepJ.z < 0) temp.prepJ = scaleproduct(-1.0,temp.prepJ);
//    cout <<FJ.IELEMENT;
//    temp.prepJ.print(cout);
//    temp.Jprep=dotproduct(temp.prepJ, FJ.J);
//    J_vec_prep_bonds.clear();


    vector<vec> vec_prep_bonds;
    vec unit_prep_bond_vec;

    vec_prep_bonds.push_back(prep_unit_vector_on_atom(FJ.IELEMENT, Atoms, Icon_atoms));
    vec_prep_bonds.push_back(prep_unit_vector_on_atom(FJ.IIELEMENT, Atoms, IIcon_atoms));
    unit_prep_bond_vec = average_vectors(vec_prep_bonds);
    temp.Jprep=dotproduct(unit_prep_bond_vec, FJ.J);
}

void calculate_desired_vectors_from_prep_unit_vec(d_flux_J &temp, vector<Atom> &Atoms, FluxJ &FJ,vector<vec> prep_unit_vecs){
    for(int i=0; i< prep_unit_vecs.size(); ++i){
        temp.J2prep.push_back(dotproduct(prep_unit_vecs[i], FJ.J));
    }
}




void read_ref_atoms_line(vector<vector<string>> &ref_atoms, string &line){
    vector<string> temp(2);
    std::istringstream iss(line);
    iss >> temp[0] >>temp[1];
    ref_atoms.push_back(temp);
}

void read_ref_bonds(vector<vector<string>> &ref_atoms, string filename){
    ifstream input;
    string line;
    input.open(filename.c_str());
    while(getline(input, line)) read_ref_atoms_line(ref_atoms, line);
}


int main(){
    vector<Atom> Atoms;
    read_geometry(Atoms, "temp.sumviz");
//    for (int i=0; i<Atoms.size(); ++i){
 //   Atoms[i].print(cout);
 //   }
    vector<FluxJ> FluxJs;
    read_extracted_FluxJ(FluxJs, "TotalJ.txt");
//    for (int i=0; i<FluxJs.size(); ++i){
//    FluxJs[i].print(cout);
//    }


/*  
  vector<d_flux_J> dfluxJ;
    for (int i = 0; i<FluxJs.size(); ++i){
        d_flux_J temp;
        temp.IELEMENT = FluxJs[i].IELEMENT;
        temp.IIELEMENT = FluxJs[i].IIELEMENT;
        vector<string> Icon_atoms,IIcon_atoms;
        Icon_atoms = connected_atoms(FluxJs[i].IELEMENT, FluxJs);
        IIcon_atoms = connected_atoms(FluxJs[i].IIELEMENT, FluxJs);
        if( Icon_atoms.size() <= 2 || IIcon_atoms.size() <= 2){
            Icon_atoms.clear();
            IIcon_atoms.clear();
            continue;
        }

        
        if( Icon_atoms.size() > 2 && IIcon_atoms.size() > 2){
            calculate_desired_vectors(temp, Atoms, FluxJs[i], Icon_atoms, IIcon_atoms);
            dfluxJ.push_back(temp);
            Icon_atoms.clear();
            IIcon_atoms.clear();

        }
        
    }
    ofstream out;
    out.open("FluxJ_prepB.txt");
*/



    vector<d_flux_J> dfluxJ;

    vector<vec> prep_unit_vecs;

    vector<vector<string>> ref_atoms;
    read_ref_bonds(ref_atoms, "ref_atoms.txt");

    prep_unit_vecs = prep_unit_vector(ref_atoms, Atoms, FluxJs);

    for (int i = 0; i<FluxJs.size(); ++i){
        d_flux_J temp;
        temp.IELEMENT = FluxJs[i].IELEMENT;
        temp.IIELEMENT = FluxJs[i].IIELEMENT;
        vector<string> Icon_atoms,IIcon_atoms;
        Icon_atoms = connected_atoms(FluxJs[i].IELEMENT, FluxJs);
        IIcon_atoms = connected_atoms(FluxJs[i].IIELEMENT, FluxJs);

//        cout<<temp.IELEMENT<< endl;
//        cout<<temp.IIELEMENT<< endl;
//        if( Icon_atoms.size() <= 2 || IIcon_atoms.size() <= 2){
//            Icon_atoms.clear();
//            IIcon_atoms.clear();
//            continue;
//        }


        
//        if( Icon_atoms.size() > 2 && IIcon_atoms.size() > 2){
            calculate_desired_vectors_from_prep_unit_vec(temp, Atoms, FluxJs[i],prep_unit_vecs);
            dfluxJ.push_back(temp);       
//        }
  //      Icon_atoms.clear();
  //      IIcon_atoms.clear();
        
    }
    ofstream out;
    out.open("FluxJ_prep2B.txt");
    out<<"\t"<<"\t";
    for (int j=0; j<ref_atoms.size(); ++j){
        out <<"\t";
        for (int k=0; k<ref_atoms[j].size();++k){
            out<<ref_atoms[j][k]<<"       ";
        }
    }
    out<<endl;

    for (int i=0; i<dfluxJ.size(); ++i){
        dfluxJ[i].print2(out);
    }
}


