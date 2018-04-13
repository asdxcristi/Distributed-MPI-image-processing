#include <iostream>
#include <mpi.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <utility>

#define INIT_PARENTS_TAG 0
#define SOBEL_TAG 1
#define MEAN_RM_TAG 2
#define STATS_AND_CLOSE_TAG 3

#define SOBEL "sobel"
#define MEAN_RM "mean_removal"

using namespace std;

struct slave_stats{
  unsigned int start_index;
  unsigned int workload_sent_size;
  unsigned int workload_recv_size;
  unsigned int sent_index;
  unsigned int recv_index;

  slave_stats():start_index(0),workload_sent_size(0),workload_recv_size(0),
  sent_index(0),recv_index(0){}
};

//loads the neighbour nodes from topology file_name
vector<int> load_neighbours(char * file_name,int rank){
  //File example format:
  /*
  0: 1 2
  1: 0 3 4 5 6
  2: 0 7 8
  3: 1
  4: 1 9 10
  5: 1
  6: 1
  7: 2
  8: 2 11
  9: 4
  10: 4
  11: 8
  */
  ifstream in_file(file_name);
  vector<int> neighs;
  int cur_node;
  char c;

  while (in_file >> cur_node >> c && c == ':'){
      string line;
      getline(in_file,line);
      if(cur_node==rank){
        istringstream is(line);
        copy(istream_iterator<int>(is),istream_iterator<int>(),
              back_inserter(neighs));
        in_file.close();
        return neighs;
    }
  }
  return neighs;
}

//load the imges list
vector<vector<string>> load_images_list(char * file_name){
  //File example format:
  /*
  1
  sobel 1.pgm 1-s.pgm
  mean_removal 1.pgm 1-m.pgm
  sobel 2.pgm 2-s.pgm
  mean_removal 2.pgm 2-m.pgm
  */
  ifstream in_file(file_name);
  vector<vector<string>> images_list;
  int nr_of_images=-1;

  in_file>>nr_of_images;

  for(int i=0;i<nr_of_images;i++){
    string fil,fi,fo;
    vector<string> cur_line;
    in_file>>fil>>fi>>fo;
    cur_line.push_back(fil);
    cur_line.push_back(fi);
    cur_line.push_back(fo);

    images_list.push_back(cur_line);
  }
  in_file.close();
  return images_list;
}

//moves the file pointer past the line containing a comment
//#define comment: string that starts with # and ends with newline
void skip_comms(ifstream& infile){
	char c;
	//save current position
	streampos oldpos = infile.tellg();
	//check if comment line
	infile>>c;
	if(c=='#'){
		while(c=='#'){
			while(c!='\n'){
				infile.get(c);
			}
			//save current position
			oldpos = infile.tellg();
			//check if next line is comment line
			//much comment galore much wow
			infile.get(c);
			if(c!='#'){
				infile.seekg(oldpos);
			}
		}
	}else{
		infile.seekg(oldpos);
	}
}

//returns a string contining a PGM file header_size
//saves the Following info:
//width -> w
//height -> h
//pixel delimitator(whitespace) -> c
string get_pgm_header(ifstream& infile,unsigned int& w,unsigned int& h,char& c){
	string inputLine = "";
	int max_val=-1;
	skip_comms(infile);
	//A "magic number" for identifying the file type. A pgm image's magic number is the two characters "P5"
	//P2 for ASCII aka in our case
  infile>>inputLine;

	skip_comms(infile);

	//A width, formatted as ASCII characters in decimal
	//A height, again in ASCII decimal
	infile>>w>>h;

	skip_comms(infile);

	//The maximum gray value (Maxval), again in ASCII decimal
	infile>>max_val;

	//Newline or other single whitespace character.
	infile.get(c);

	int header_size = infile.tellg();
	int dummy;
	infile>>dummy;
	infile.get(c);//c is now pixel delimitator

	//Save the header
	char * header_buffer=new char[header_size+1];
  header_buffer[header_size]='\0';
	infile.seekg(0,infile.beg);
	infile.read(header_buffer, header_size);

	string header=header_buffer;

  delete[] header_buffer;
	return header;
}

//returns the tag coresponding to the filter
//check #define section or README for numbers
int get_filter_tag(string filter){
  int tag=-1;

  if(filter==SOBEL){
    tag=SOBEL_TAG;
  }else if(filter==MEAN_RM){
    tag=MEAN_RM_TAG;
  }

  return tag;
}

//applies filter to the matrix of pixels
//check #define section or README for numbers
//supported filters: SOBEL; MEAN_RN
vector<int> solve_filter(const vector<int>& raw,const int& width, const int& height,const int& filter){
  vector<int> result;
  result.resize((width-2)*(height-2));
  //kernel positions relative to center
  vector<pair<int,int> > kernel={
  {-1,-1},{-1,0},{-1,1},
  {0,-1},{0,0},{0,1},
  {1,-1},{1,0},{1,1}};

  //sobel filter
  vector<int> sobel_kernel={
  1,0,-1,
  2,0,-2,
  1,0,-1};

  //mean removal filter
  vector<int> mean_rm_kernel={
    -1,-1,-1,
    -1,9,-1,
    -1,-1,-1};

    for(int i=0;i<height-2;i++){
     for(int j=0;j<width-2;j++){
       switch (filter) {
         case SOBEL_TAG:
         for(unsigned int it=0;it<kernel.size();it++){
           result[i*(width-2)+j]+=raw[(i+1+kernel[it].first)*width+(j+1+kernel[it].second)]*sobel_kernel[it];
         }
          result[i*(width-2)+j]/=1;
          result[i*(width-2)+j]+=127;
          break;

         case MEAN_RM_TAG:
         for(unsigned int it=0;it<kernel.size();it++){
           result[i*(width-2)+j]+=raw[(i+1+kernel[it].first)*width+(j+1+kernel[it].second)]*mean_rm_kernel[it];
         }

         result[i*(width-2)+j]/=1;
         result[i*(width-2)+j]+=0;
         break;

         default:
          result[i*(width-2)+j]=raw[(i+1)*width+(j+1)];
          break;
       }

       if(result[i*(width-2)+j]>255){
         result[i*(width-2)+j]=255;
       }
       if(result[i*(width-2)+j]<0){
        result[i*(width-2)+j]=0;
      }
    }//for rows
  }//for cols
  return result;
}


int main(int argc,char *argv[]){
  //safeguard
  if(argc!=4){
    cout<<"[ERROR]Usage : mpirun -np N ./filtru topologie.in imagini.in statistica.out \n";
    return -1;
  }

  string stats_file_name=argv[3];


	MPI_Init(&argc, &argv);
  int rank,nProcesses;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

  vector<int> neighs;

  neighs=load_neighbours(argv[1],rank);

  int parent=-1;
  int probe=1337;
  int tag=1337;
  //Send probe to init parents
  if(rank==0){
    tag=INIT_PARENTS_TAG;
    for(auto i:neighs){
      MPI_Send(&probe, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
    }

    for(auto i:neighs){
      MPI_Recv(&probe, 1, MPI_INT, i, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    //Start processing images
    vector<vector<string>> images_list=load_images_list(argv[2]);
    for(auto cur_image_stuff:images_list){
      string input_file,output_file,filter;

      filter=cur_image_stuff[0];
      input_file=cur_image_stuff[1];
      output_file=cur_image_stuff[2];

      auto filter_tag=get_filter_tag(filter);

      unsigned int w = 0, h = 0;
    	char c;
      ifstream infile(input_file);

    	string header=get_pgm_header(infile,w,h,c);

      unsigned int borded_width=w+2;
      unsigned int borded_height=h+2;

      vector<unsigned int> array;
      array.resize(borded_height*borded_width);

    	for(unsigned int i=0;i<borded_height;i++)
    		for(unsigned int j=0;j<borded_width;j++)
    			array[i*borded_width+j]=0;

      //following lines : data
      for(unsigned int row=1;row<borded_height-1;row++)
        for(unsigned int col=1;col<borded_width-1;col++)
    			infile>>array[row*borded_width+col];

    	infile.close();

      slave_stats stats[neighs.size()];
      //compute the workload 4 every child
      int workload_size=(borded_height-2)/neighs.size();
      int send_from_index=0;
      int recv_to_index=0;

      //send work to all childs
      for(unsigned int i=0;i<neighs.size();i++){
        MPI_Send(&borded_width, 1, MPI_INT,neighs[i], filter_tag, MPI_COMM_WORLD);

        if(i==neighs.size()-1 && workload_size*neighs.size() != borded_height-2){
          workload_size+=borded_height-2-workload_size*neighs.size();
        }

        stats[i].sent_index=send_from_index;
        stats[i].recv_index=recv_to_index;
        stats[i].workload_sent_size=borded_width*(workload_size+2);
        stats[i].workload_recv_size=workload_size*w;
        send_from_index+=borded_width*workload_size;
        recv_to_index+=workload_size*w;

        MPI_Send(&array[stats[i].sent_index],stats[i].workload_sent_size,MPI_INT,neighs[i],filter_tag,MPI_COMM_WORLD);
      }

      array.clear();
      array.resize(w*h);
      //receive proccesed pixels
      for(unsigned int i=0;i<neighs.size();i++){
        MPI_Recv(&array[stats[i].recv_index], stats[i].workload_recv_size, MPI_INT, neighs[i], filter_tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }

      //write resulted image to file
      ofstream of(output_file);
    	of<<header;

  	 for(unsigned int row=0;row<h;row++)
       for(unsigned int col=0;col<w;col++)
  			of<<array[row*w+col]<<c;

    	of.close();
      array.clear();
    }//for every images


    vector<int> stats;
    stats.resize(nProcesses);
    vector<int> recvd_stats;
    recvd_stats.resize(nProcesses);

    //Send exit signal
    for(auto neigh:neighs){
      MPI_Send(&stats[0], stats.size(), MPI_INT,neigh, STATS_AND_CLOSE_TAG, MPI_COMM_WORLD);
    }
    //mergs all statistics
    for(auto neigh:neighs){
      MPI_Recv(&recvd_stats[0],recvd_stats.size(),MPI_INT, neigh, STATS_AND_CLOSE_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      for(unsigned int i=0;i<stats.size();i++){
        stats[i]|=recvd_stats[i];
      }
    }

    //Write stats 2 file
    ofstream stats_file(stats_file_name);
    for(unsigned int i=0;i<stats.size();i++){
      stats_file<<i<<": "<<stats[i]<<endl;
    }
    stats_file.close();
  }//rank 0

  if(rank!=0){
    MPI_Status status;
    int size;
    unsigned int width;
    unsigned int nr_of_processed_lines=0;
    vector<int> stats;
    vector<int> recvd_stats;

    while(true){
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      //check message type using tag
      switch (status.MPI_TAG) {
        case INIT_PARENTS_TAG:
          MPI_Recv(&probe, nProcesses, MPI_INT, MPI_ANY_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,&status);

          parent=status.MPI_SOURCE;
          //aruncam la gunoi parintele din lista de vecini
          neighs.erase(remove(neighs.begin(),neighs.end(),parent),neighs.end());

          //forward probe to neighs
          for(int i:neighs){
            MPI_Send(&probe, 1, MPI_INT, i, status.MPI_TAG, MPI_COMM_WORLD);
          }
          //await resp from neigh
          for(int i:neighs){
            MPI_Recv(&probe, 1, MPI_INT, i, status.MPI_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          }
          MPI_Send(&probe, 1, MPI_INT,parent, status.MPI_TAG, MPI_COMM_WORLD);
          break;

        case STATS_AND_CLOSE_TAG:
          MPI_Get_count(&status, MPI_INT, &size);
          stats.resize(size);
          recvd_stats.resize(size);

          MPI_Recv(&recvd_stats[0], recvd_stats.size(), MPI_INT, parent, status.MPI_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

          //if intermediary node, announce others
          if(neighs.size()>0){
            for(auto neigh:neighs){
              MPI_Send(&recvd_stats[0], recvd_stats.size(), MPI_INT, neigh, status.MPI_TAG, MPI_COMM_WORLD);
            }
            //recv all statistics
            for(auto neigh:neighs){
              MPI_Recv(&recvd_stats[0], recvd_stats.size(), MPI_INT, neigh, status.MPI_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

              //mergs all statistics
              for(unsigned int i=0;i<stats.size();i++){
                stats[i]|=recvd_stats[i];
              }
            }
          }

          stats[rank]=nr_of_processed_lines;

          MPI_Send(&stats[0], stats.size(), MPI_INT,parent, status.MPI_TAG, MPI_COMM_WORLD);

          MPI_Finalize();
          return 0;

        case MEAN_RM_TAG://HUE HUE HUE
        case SOBEL_TAG:
          MPI_Get_count(&status, MPI_INT, &size);
          if(size==1){
            MPI_Recv(&width, size, MPI_INT, parent, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }else{
            unsigned int height=size/width;
            vector<int> chunk;
            chunk.resize(size);

            for(unsigned int i=0;i<height;i++)
              for(unsigned int j=0;j<width;j++)
                chunk[i*width+j]=0;

            MPI_Recv(&chunk[0], chunk.size(), MPI_INT, parent, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if(neighs.size()>0){
              slave_stats slave_info[neighs.size()];

              unsigned int workload_size=(height-2)/neighs.size();
              unsigned int send_from_index=0;
              unsigned int recv_to_index=0;

              for(unsigned int i=0;i<neighs.size();i++){
                MPI_Send(&width, 1, MPI_INT,neighs[i], status.MPI_TAG, MPI_COMM_WORLD);

                if(i==neighs.size()-1 && workload_size*neighs.size() < height-2){
                  workload_size+=height-2-workload_size*neighs.size();
                }

                slave_info[i].sent_index=send_from_index;
                slave_info[i].recv_index=recv_to_index;
                slave_info[i].workload_sent_size=width*(workload_size+2);
                slave_info[i].workload_recv_size=workload_size*(width-2);
                send_from_index+=width*workload_size;
                recv_to_index+=slave_info[i].workload_recv_size;


                MPI_Send(&chunk[slave_info[i].sent_index],slave_info[i].workload_sent_size,MPI_INT,neighs[i],status.MPI_TAG,MPI_COMM_WORLD);
              }

              chunk.clear();
              chunk.resize((height-2)*(width-2));

              for (unsigned i = 0; i < neighs.size(); i++) {
                MPI_Recv(&chunk[slave_info[i].recv_index], slave_info[i].workload_recv_size, MPI_INT, neighs[i], status.MPI_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
              }

              MPI_Send(&chunk[0],chunk.size(),MPI_INT,parent,status.MPI_TAG,MPI_COMM_WORLD);

            }else{//if end
              //leaf working..poor child
              vector<int> res=solve_filter(chunk,width,height,status.MPI_TAG);
              nr_of_processed_lines+=height-2;

              MPI_Send(&res[0],res.size(),MPI_INT,parent,status.MPI_TAG,MPI_COMM_WORLD);

              res.clear();
            }
            chunk.clear();
          }
          break;
      }//switch tag

    }//while

  }


  MPI_Finalize();
  return 0;
}
