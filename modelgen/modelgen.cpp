
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <sstream>
using namespace std;

unordered_map<string, bool> token_map = {
	{":", 1},
	{";", 1},
	{"reversible", 1},
	{"Main", 1},
	{"Modules", 1},
	{"DPD", 1},
	{"AdditionalRules", 1},

};

unordered_map<char,bool> break_map = {
	{' ', 1},
	{'.', 1},
	{',', 1},
	{';', 1},
	{':', 1},
	{'(', 1},
	{')', 1},
	{'\t', 1},
	{'{', 1},
	{'}', 1},
	{'+', 1},
	{'!', 1},
	{'=', 1}
};


vector<string> tokenize() {
	
	ifstream cstar_file("../model/model.cstar");
	string line;
	vector<string> tokens = {};
	string buffer = "";
	bool comment_flag = false;
	bool str_flag = false;
		
	while(getline(cstar_file, line)) {
		comment_flag = false;
		for(int i = 0; i < line.size(); i++) {
			if(line[i] == '\"') {
				str_flag = !str_flag;
				if(buffer.size()) tokens.push_back(buffer);
				buffer = "";
				continue;
			}
			
			if(str_flag) {
				buffer += line[i];
				continue; 
			}

		
			
			if(line[i] == '#') comment_flag = true;
			if(comment_flag) continue;

			if(break_map.find(line[i]) != break_map.end()) {
				// deal with single tokens
				if(buffer == "" and line[i] != ' ' and line[i] != '\t') {
					buffer += line[i];	
					if(buffer.size()) tokens.push_back(buffer);
					buffer = "";
					continue;

				}
				
				if(buffer != " " and buffer.size()) tokens.push_back(buffer);
				//cout << buffer << ":" << buffer.size() << endl;
				buffer = "";
				

				if(line[i] != ' ' and line[i] != '\t') tokens.push_back(buffer+line[i]);
				continue;
			}

			buffer += line[i];
		}
			
	}

	return tokens;
}


// parsing functions
vector<vector<string>> parse_block(vector<string> tokens) {
	vector<vector<string>> block_matrix = {};
	vector<string> buffer_vec = {};
	for(auto& token: tokens) {
		if(token == ";") {
			block_matrix.push_back(buffer_vec);
			buffer_vec.clear();
			continue;
		}

		buffer_vec.push_back(token);
		
	}

	return block_matrix;

}


vector<pair<string, vector<string>>> parse_program(vector<string> tokens) {
	vector<pair<string, vector<string>>> program_matrix = {}; 
	vector<string> buffer_vec = {};
	string current_type = "";
	int brace_counter = 0;
	for(auto& token: tokens) {
		if(token == "{") {
			if(!brace_counter) {
				current_type = buffer_vec[buffer_vec.size()-1];
				buffer_vec.clear();
			}
	
			brace_counter++;
			continue;
		
		} 

		if(token == "}") {
			brace_counter--;
			if(!brace_counter) {
				pair<string, vector<string>> subprogram(current_type, buffer_vec);
				program_matrix.push_back(subprogram);
				buffer_vec.clear();
				continue;
			}
		}

		buffer_vec.push_back(token);	
	}
	
	return program_matrix;

} 


struct Module {
	string name;
	string inhibitor;
	string ligand;
	string activation_site;
};


int get_index(vector<string> vec, string str) {
	auto it = find(vec.begin(), vec.end(), str);
	if(it != vec.end()) {
		return (int)(it - vec.begin());
	}

	return -1;
}



vector<Module> get_modules(vector<string> module_block) {
	// get modules
	vector<vector<string>> modules = parse_block(module_block);
	vector<Module> out_modules = {};
	for(auto& mod : modules) {
		Module new_mod;
		new_mod.name = mod[0];
		//cout << new_mod.name << endl;
		int inhibit_index = get_index(mod, ":") + 1;
		new_mod.inhibitor = inhibit_index == 0 ? "" :  mod[inhibit_index];
		
		int ligand_index = get_index(mod, ",") + 1;
		new_mod.ligand = ligand_index == 0 ?  "" : mod[ligand_index];

		int active_index = get_index(mod, "->") + 1;
		new_mod.activation_site = active_index == 0 ? "" :  mod[active_index];
		out_modules.push_back(new_mod);

	}

	return out_modules;


}


class Model {
	public:
	
		Model() {
			network.open("/Users/atamb/compbio/consensus_MCF10A.csv");
			file.open("model2test.bngl");	
			file << "begin model\n";
			net_modules = {};
		}

		~Model() {
			file.close();
		}

		void add_dpd(string dpd) {
			this->dpd_vals.push_back(dpd);
		}


		void add_modules(vector<Module> core_mods) {
			this->net_modules = core_mods;
			file << "begin molecule types\n";
			for(Module mod : core_mods) {
				string lig = mod.ligand != "" ? "ligand," : "";
				string inh = mod.inhibitor != "" ? "inhibitor," : "";
				string as = mod.activation_site != "" ? mod.activation_site + "~0~p" : "";
				file << mod.name << "(" << lig << inh << as << ")\n";
				
				if(inh != "") file << mod.inhibitor << "(" << mod.name << "Rec" << ")\n";

			}

			file << "end molecule types\n";
		}


		void generate_dpd_force() {
						


		}

		void gen_observables() {
				

		}


		void gen_funcs() {
			string line;
			int counter = 0;
			vector<vector<string>> matrix = {};
			while(getline(network, line)) {
				istringstream ss(line);
				string token;
				vector<string> row = {};
				vector<string> module_order = {};
				while(std::getline(ss, token, ',')) {
				    if(token == "" or token == " ") continue;
				    //cout << token << endl;
				    if(counter == 0) module_order.push_back(token);
				    else row.push_back(token);
				
				}

				if(counter != 0) matrix.push_back(row);
				
				counter = 1;

			}


			for(int i = 0; i < net_modules.size() + 2; i++) {
				for(int j = 0; j < net_modules.size() + 2; j++) {
					//cout << i << j << endl;
					cout << matrix[i][j];

				}

				cout << endl;
				
			}
			
		}

		
	private:	
		ifstream network;
		ofstream file;
		vector<Module> net_modules;
		vector<string> dpd_vals;
};



int main() {
	Model new_net_model;

	vector<string> tokens_out = tokenize();	
	for(int i = 0; i < tokens_out.size(); i++) {
	//	cout << tokens_out[i] << " , ";

	}

	vector<pair<string, vector<string>>> parsed_body = parse_program(tokens_out);
	for(auto& part : parsed_body) {
		if(part.first == "modules") {
			vector<Module> model_mods = get_modules(part.second);
			new_net_model.add_modules(model_mods);

		}
	}

	new_net_model.gen_funcs();
	


}
