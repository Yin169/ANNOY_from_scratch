#include <iostream>
#include <vector>
#include <cmath>

struct Node{
    std::vector<float> val;
    int side;
};

typedef std::pair<Node, Node> PII;

class AnnoyIndex{
    public:

    AnnoyIndex* left;
    AnnoyIndex* right;
    std::vector<float> mid_point, pq_line;
    std::vector<Node> leaf_node;

    AnnoyIndex(AnnoyIndex* _left, AnnoyIndex* _right, int _dim, int _datasize, std::vector<Node> _leaf_node) : left(_left), right(_right), dim(_dim), datasize(_datasize), leaf_node(_leaf_node) {
        for (size_t i=0; i<dim; i++){
            mid_point.push_back(0.0);
            pq_line.push_back(0.0);
        }
    };
    ~AnnoyIndex() {};

    AnnoyIndex* fit(std::vector<Node> nodes);
    float get_distance(Node a, Node b);
    float dot_product(Node a, Node b);
    PII TwoMean(std::vector<Node> nodes);
    Node Normalize(Node a);
    
    protected:
    const static float eps = 1e-6;
    const static size_t num_iter = 100;
    const static size_t max_ind = 12;
    size_t dim;
    size_t datasize;
};

class AnnoyInterface{
    public:
    AnnoyIndex* root;

    AnnoyInterface(std::vector<Node> nodes, int _dim) : dataset(nodes), dim(_dim) {};
    ~AnnoyInterface() {};
    bool build();
    void get_nns_vector(Node query, AnnoyIndex* root);
    
    protected:
    std::vector<Node> dataset;
    size_t dim; 
};

Node AnnoyIndex :: Normalize(Node a){
    float res = 0;
    for (size_t i=0; i<dim; i++){
        res += a.val[i]*a.val[i];
    }
    res = std::sqrt(res);
    for (size_t i=0; i<dim; i++){
        a.val[i] /= res;
    }
    return a;
}

float AnnoyIndex :: get_distance(Node a, Node b){
    float res = 0;
    for (size_t i=0; i<dim; i++){
        res += (a.val[i]-b.val[i])*(a.val[i]-b.val[i]);
    }
    return std::sqrt(res);
}

float AnnoyIndex :: dot_product(Node a, Node b){
    // a = Normalize(a), b = Normalize(b);
    float res = 0;
    for (size_t i=0; i<dim; i++){
        res += a.val[i]*b.val[i];
    }
    return res;
}

void showNode(Node a, size_t dim){
    for (size_t i=0; i<dim; i++){
        std::cout << a.val[i] << " ";
    }
    std::cout << std::endl;
}

PII AnnoyIndex :: TwoMean(std::vector<Node> nodes){
    std::vector<float> center1, center2;
    for (size_t i=0; i<dim; i++){
        center1.push_back((std::rand()%nodes.size()/(float)nodes.size()) + 0.5);
        center2.push_back((std::rand()%nodes.size()/(float)nodes.size()));
    }
    Node cen1 = {center1, -1};
    Node cen2 = {center2, -1};
    // showNode(cen1, dim);
    // showNode(cen2, dim);
    size_t max_iter_number = num_iter;
    while (max_iter_number--){
        int cnt1 = 0, cnt2=0;
        Node prev_cen1 = cen1;
        Node prev_cen2 = cen2;
        std::fill(cen1.val.begin(), cen1.val.end(), 0.0);
        std::fill(cen2.val.begin(), cen2.val.end(), 0.0);
        for (Node node : nodes){
            float dist1 = get_distance(prev_cen1, node);
            float dist2 = get_distance(prev_cen2, node);
            for (size_t i=0; i<dim; i++){
                if (dist1<dist2) {
                    cen1.val[i]+=node.val[i];
                    cnt1++;
                }
                else {
                    cen2.val[i]+=node.val[i];
                    cnt2++;
                }
            }
        }

        for (size_t i=0; i<dim; i++){
            cen1.val[i] /= (cnt1+1e-6);
            cen2.val[i] /= (cnt2+1e-6);
        }
        // showNode(cen1, dim);
        // showNode(cen2, dim);
    }
    PII res(cen1, cen2);
    return res;
}

AnnoyIndex* AnnoyIndex :: fit(std::vector<Node> nodes){
    std::cout << nodes.size() << std::endl;
    if (nodes.size() <= max_ind) return new AnnoyIndex(left=nullptr, right=nullptr, dim=dim, datasize=nodes.size(), leaf_node=nodes);
    PII cens = TwoMean(nodes);
    Node p = cens.first, q = cens.second;
    float a = dot_product(p, q);
    // std::cout << "pq " << a << std::endl;
    if (a < 0) std::swap(p, q); 
    // showNode(p, dim);
    // showNode(q, dim);
    for (size_t i=0; i<dim; i++){
        mid_point[i] = (p.val[i]+q.val[i])/2;
        pq_line[i] = (p.val[i]-q.val[i]);
    }
    Node pq = {pq_line, -1};
    std::vector<Node> nodes1;
    std::vector<Node> nodes0;
    for (Node node : nodes){
        std::vector<float> tmp;
        for (size_t i=0; i<dim; i++){
            tmp.push_back(node.val[i]-mid_point[i]);
        }
        Node tmp_node = {tmp, -1};
        float a = dot_product(tmp_node, pq);
        if (a > 0 && std::rand()%2==1) {
            node.side = 1;
            nodes1.push_back(node);
        }
        else {
            node.side = 0;
            nodes0.push_back(node);
        }
    }
    std::vector<Node> empty;
    return new AnnoyIndex(left=fit(nodes1), right=fit(nodes0), dim=dim, datasize=nodes.size(), leaf_node=empty);
}


bool AnnoyInterface :: build (){
    std::vector<Node> empty;
    AnnoyIndex annoy(nullptr, nullptr, dim, dataset.size(), empty);
    root = annoy.fit(dataset);
    return root != nullptr ? true : false;
}

void AnnoyInterface :: get_nns_vector(Node query, AnnoyIndex* root){
    if (root->leaf_node.size()>0){
        for (Node node : root->leaf_node) {
            node = root->Normalize(node), query = root->Normalize(query);
            showNode(node, dim);
            std::cout << root->dot_product(node, query) << std::endl;
        }
        return;
    }
    Node pq = {root->pq_line, -1};
    
    std::vector<float> tmp;
    for (size_t i=0; i<dim; i++){
        tmp.push_back(query.val[i] - root->mid_point[i]);
    }
    Node tmp_node = {tmp, -1};
    float a = root->dot_product(tmp_node, pq);
    if (a > 0 ) get_nns_vector(query, root->left);
    else get_nns_vector(query, root->right);
}
