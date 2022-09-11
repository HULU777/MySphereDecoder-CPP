#include "mwvc.h"
#include <sstream>
#include <fstream>

int main(int argc, char *argv[])
{
    uint seed;
    
    if (argc == 1)
    {
        cout << "FastWVC - a Minimum Weighted Vertex Cover Problem solver." << endl;
        cout << "Usage: ./mwvc [Graph file] [Seed] [Cutoff time] [CC mode]" << endl;
        return 1;
    }

    if (argc < 5)
    {
        cerr << "Missing argument(s)." << endl;
        cout << "Usage: ./mwvc [Graph file] [Seed] [Cutoff time] [CC mode]" << endl;
        return 1;
    }

    stringstream ss;  // char to stringstream to int
    ss << argv[2];
    ss >> seed;
    ss.clear();
    ss << argv[3];
    ss >> cutoff_time;
    ss.clear();
    ss << argv[4];
    ss >> mode;
    ss.clear();

    // if (BuildInstance(argv[1]) != 0)
    // {
    //     cerr << "Open instance file failed." << endl;
    //     return 1;
    // }

    if (seed < 0U || seed > ~0U)
    {
        seed = 10;
    }

    if (cutoff_time < 0 || cutoff_time > (int)(~0U>>1))
    {
        cutoff_time = 1000;
    }

    if (mode < 0 || mode > 3)
    {
        mode = 0;
    }

    srand(seed);

    cout << argv[1];

    start = chrono::steady_clock::now();

    ConstructVC(); 
    LocalSearch();
    int v;
    if (CheckSolution() == 1)
    {
        cout << ", " << best_weight << ", " << best_comp_time << endl;

        ofstream file;
        file.open("/home/hliucp/Desktop/VertexCover/MySphereDecoder/okk.txt");  //open a file  VCforSD/
        cout << v_num << endl;
        for (v = 1; v < v_num+1; v++)
        {
            // cout << "," << best_v_in_c[v] << endl;
            file<< *(best_v_in_c+v) << endl;       //write to it
            // cout<< *(best_v_in_c+v)<< endl;       //result = Hello file
        }
        file.close();           //close it
    }
    else
    {
        cout << ", the solution is wrong." << endl;
    }

    FreeMemory();

    return 0;
}
