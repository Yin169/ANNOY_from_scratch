#include <iostream>
#include <vector>
#include "src/annoy.hh"

int main(){
    std::cout << "Starting..." << std::endl;
    std::vector<Node> example;
    size_t dim = 2;
    size_t number = 1000;
    for (int i =0; i<number; i++){
        std::vector<float> tmp, mtmp;
        for (int j=0; j<dim; j++){
            float val = (std::rand()%number) /(float)number;
            tmp.push_back(val);
            mtmp.push_back(-1.0*val);
        }
        Node node = {tmp, -1};
        Node mnode = {mtmp, -1};
        example.push_back(node);
        example.push_back(mnode);
    }
    std::cout << "Create interface..." << std::endl;
    AnnoyInterface annoy_interface(example, dim);
    std::cout << annoy_interface.build() << std::endl;

    std::cout << "ANN Searching..." << std::endl;
    showNode(example[0], dim);
    std::cout << "........" << std::endl;
    annoy_interface.get_nns_vector(example[0], annoy_interface.root);
    return 0;
}