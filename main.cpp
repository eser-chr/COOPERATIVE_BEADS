#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>  
#include <random>
#include <chrono>
#include <algorithm> 
// #include <omp.h>

namespace program_options {

struct Options {
    size_t X; 
    size_t T; 
    size_t N; 
    size_t SIMS; 
    int k_on; // We use divider to adjust prob to bind
    int k_off; // Similarly for the unbind
    int step_prob; // 1/10 
    const int divider=1000000;
    size_t blocks;
    int period;
    // int period=T/blocks; // Number of blocks. At each block we will write the result
    
    void print() const {
        std::cout<<"############################### \n";
        std::cout<<"Num of Sims "<<SIMS<<"\n";
        std::cout<<"Num of sides "<<X<<"\n";
        std::cout<<"Num of timesteps "<<T<<"\n";
        std::cout<<"Num of kinesins "<< N<<"\n";
        std::cout<<"Num of blocks "<< blocks<<"\n";
        std::cout<<"binding prob "<<static_cast<double>(k_on)/static_cast<double>(divider)<<"\n";
        std::cout<<"unbinding prob "<<static_cast<double>(k_off)/static_cast<double>(divider)<<"\n";
        std::cout<<"step prob "<<static_cast<double>(step_prob)/static_cast<double>(divider)<<"\n";
        std::cout<<"###############################"<<std::endl;
    }
};

auto parse(int argc, char *argv[]) {
    if (argc != 9)
        throw std::runtime_error("unexpected number of arguments");
    Options opts;
    
    if (std::sscanf(argv[1], "%zu", &opts.X) != 1)
        throw std::runtime_error("invalid parameter for (X) number of sides");
    if (std::sscanf(argv[2], "%zu", &opts.T) != 1)
        throw std::runtime_error("invalid parameter for (T) number of steps");
    if (std::sscanf(argv[3], "%zu", &opts.N) != 1)
        throw std::runtime_error("invalid parameter for (N) number of kinesins");
    if (std::sscanf(argv[4], "%zu", &opts.SIMS) != 1)
        throw std::runtime_error("invalid parameter for (SIMS) number of simulations");


    if (std::sscanf(argv[5], "%d", &opts.k_on) != 1)
        throw std::runtime_error("invalid parameter for (k_on) probabitlity to bind");
    if (std::sscanf(argv[6], "%d", &opts.k_off) != 1)
        throw std::runtime_error("invalid parameter for (k_off) probability to unbind");
    if (std::sscanf(argv[7], "%d", &opts.step_prob) != 1)
        throw std::runtime_error("invalid parameter for (step_prob) probability to step forward");

    if (std::sscanf(argv[8], "%zu", &opts.blocks) != 1)
        throw std::runtime_error("invalid parameter for (blocks) number of blocks");
    opts.period=opts.T/opts.blocks;

    return opts;
}

} // namespace program_options



struct Kinesin{
public:
    int position;
    bool isBound;
    Kinesin()=default;
    Kinesin(int position, bool isBound): position(position), isBound(isBound){};
};



class Simulation{
    int id;
    int seed;
    program_options::Options opts;
    std::mt19937 engine;
    std::vector<int> DATA;
    std::vector<Kinesin> kinesin_list;
    std::vector<bool> grid;
    // std::uniform_int_distribution<int> placer; //(0,opts.X);
    std::uniform_int_distribution<int> roller;//(0,opts.divider);

    struct Initialization{
        static void set_to_zero(Simulation* sim){    // Cleaning any previous data from other simulations.
            sim->kinesin_list.clear();
            std::fill(sim->grid.begin(), sim->grid.end(), false);	
        }

        static void initialize(Simulation* sim){
            for(int i=0; i<sim->opts.N; i++){		
                sim->kinesin_list[i] = Kinesin(i+3000, true);
                sim->grid[i+3000]=true;
            }	
        }
    };

    class Move{
        Simulation* sim;
        int kin;
        int position;
        bool type;
        bool grid_x;
        bool grid_x_;
    
    public:
        Move(Simulation* sim, int& kin):sim(sim), kin(kin){
            position = sim->kinesin_list[kin].position;
            type = sim->kinesin_list[kin].isBound;
            grid_x = sim->grid[position];
            grid_x_ = sim->grid[position+1];
        }

        void set_bind(){
            int x= sim->kinesin_list[kin].position;
            sim->kinesin_list[kin].isBound=true;
            sim->grid[x]=true;	
        }

        void set_unbind(){ //set kinesin to unbind
            int x= sim->kinesin_list[kin].position;
            sim->kinesin_list[kin].isBound=false;
            sim->grid[x]=false;
        }

        void set_move(){ //set kinesin to unbind
            int x= sim->kinesin_list[kin].position;
            sim->grid[x]=false;
            sim->kinesin_list[kin].position=x+1;		
            sim->grid[x+1]=true;
        }


        int prob_to_bind(){
            int x = sim->kinesin_list[kin].position;
            bool type =sim->kinesin_list[kin].isBound;

            if (type==1 || sim->grid[x]==false){ // If its bound or another kin is bound at that place
                return 0;
            }else{
                return sim->opts.k_on;
            }
        }

        int prob_to_unbind(){
            int x = sim->kinesin_list[kin].position;
            bool type =sim->kinesin_list[kin].isBound;

            if (type==0){
                return 0;
            }else{
                return sim->opts.k_off;
            }
        }

        int prob_to_step(){
            int x = sim->kinesin_list[kin].position;
            bool type =sim->kinesin_list[kin].isBound;

            if (type==0 || sim->grid[x+1]==true || x+1>=sim->opts.X ){
                return 0;
            }else{
                return sim->opts.step_prob;
            }

        }

    };

    

    void Gilespie(int kin, std::mt19937& engine){
        int random_number =roller(engine);
        Move move(this, kin);
        int binding = move.prob_to_bind();
        int unbinding = binding + move.prob_to_unbind();
        int step = unbinding + move.prob_to_step();

        if(random_number<binding){
            move.set_bind();
        }else if (random_number<unbinding){
            move.set_unbind();
        }else if(random_number<step){
            move.set_move();
        }

    }



public:
    Simulation(program_options::Options opts, int seed, int id):opts(opts), seed(seed), id(id){
        DATA.reserve(opts.N*opts.blocks);
        kinesin_list.reserve(opts.N);
        grid.reserve(opts.X);
        engine.seed(seed);
        roller.param(std::uniform_int_distribution<int>::param_type(0, opts.divider));
    }

    void save_step(int& step){
    if(step%opts.period==0){
			int i=step/opts.period;			
			for(int Kin=0; Kin<opts.N;Kin++){	
                DATA[i*opts.N+Kin]=kinesin_list[Kin].position;				
			}
        }	
    }

    void write_to_file(){ 
        std::string name = "./data"+std::to_string(id)+".csv";
        std::ofstream file;
        file.open(name, std::ofstream::trunc);
        for(int i =0; i<opts.blocks-1; ++i){
                file<<i<<",";
            }
        file<<opts.blocks-1<<std::endl;
        for(int kin = 0; kin<opts.N; ++kin){
            for(int i =0; i<opts.blocks -1; ++i){
                file<<DATA[i*opts.N+kin]<<",";
            }
            file<<DATA[(opts.blocks-1)*opts.N + kin]<<std::endl;
        }
    }
    

    void run(){
        Initialization::initialize(this);
        std::vector<int> vec(opts.N);
        std::iota(vec.begin(), vec.end(), 1);
        for(int step=0; step<opts.T; step++){
            save_step(step);
            std::shuffle(vec.begin(), vec.end(),engine);
            for(int kin:vec){
                Gilespie(kin-1,engine);			
            }	
        }
        write_to_file();
    }

};





int main(int argc, char *argv[]){
    std::cout << std::endl;
  // parse args
    auto opts = program_options::parse(argc, argv);
    opts.print();
    std::uniform_int_distribution<int> seeder(0, 1000);
    std::mt19937 engine(10);



#ifdef CALC_TIME
    auto start_time = std::chrono::high_resolution_clock::now();
#endif
	// unsigned int seed =15;	

	#pragma omp parallel for
	for (int i=0; i<opts.SIMS;i++){
        
		Simulation sim(opts, seeder(engine), i);
        sim.run();
	}


#ifdef CALC_TIME
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Runtime: " << duration.count() << " milliseconds" << std::endl;
#endif
}
































