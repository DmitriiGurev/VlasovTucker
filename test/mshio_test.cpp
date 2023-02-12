#include <mshio/mshio.h>

#include <iostream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    string fileName = argv[1];

    mshio::MshSpec spec = mshio::load_msh(fileName);

    cout << "Nodes:\n";

    cout << right << fixed << setprecision(4);
    for (int i = 0; i < spec.nodes.num_nodes; i++)
	{
		cout << setw(5) << spec.nodes.entity_blocks[0].tags[i] << " "
			 << setw(7) << spec.nodes.entity_blocks[0].data[3 * i] << " "
			 << setw(7) << spec.nodes.entity_blocks[0].data[3 * i + 1] << " "
			 << setw(7) << spec.nodes.entity_blocks[0].data[3 * i + 2] << "\n";
	}
}