#include <iostream>
#include <vector>
#include "src/annoy.hh"

int main(){
    std::cout << "Starting..." << std::endl;
    std::vector<Node> example;
    size_t dim = 10;
    size_t number = 1000;
    for (int i =0; i<number; i++){
        std::vector<float> tmp;
        for (int j=0; j<dim; j++){
            tmp.push_back(std::rand()%number /(float)number);
        }
        Node node = {tmp, -1};
        example.push_back(node);
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