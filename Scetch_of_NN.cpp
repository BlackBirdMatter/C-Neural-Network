#include <iostream>
#include <cmath>
#include <algorithm>
#include<array>
using namespace std;

class CrossEntropy{
    private:
        double probs[512]; 
        double real[512]; 
        double arrSize;

    public:
        void setParams(double probs[], double real[], int arrSize){
            for (int i = 0; i < arrSize; i++){
                this->probs[i] = probs[i];
                this->real[i]  = real[i];
            }
            this-> arrSize  =  double(arrSize);
        }
        double fit(){
            double summ = 0;
            for (int i = 0; i < arrSize; i++){
                summ += (real[i] * log(probs[i])) + ((1 - real[i]) * log((1 - probs[i])));
                
            }
            double crossentropy = -((1 / arrSize) * summ);
            return crossentropy;
        }
};

class NeuralNetwork{
    private:
        int input_dim;
        int output_dim;
    

};

class Neuron{
    private:
        int input_size;
        
        
    public:
        double last_output;
        double weights[1024];
        void setParams(double weights[], double intercept, int input_size){
            for (int i = 0; i < input_size; i++){
                this->weights[i] = weights[i];
            }
            this->weights[input_size+1]  =  intercept;
            this-> input_size = input_size;
        };
        double fit(double input_values[]){
            double summ = 0.0;
            for (int i = 0; i < input_size; i++){
                summ += weights[i] * input_values[i];
            }
            summ += weights[input_size+1];
            this->last_output = summ;
            return summ;
        }

};

class Sigmoid{
    public: 
        double fit(double z){
            return 1 / (1 + exp(-z));
        };
};

class Tanh{
    public: 
        double fit(double z){
            return (exp(z) - exp(-z))/(exp(z) + exp(-z));
        };
};

class Softmax{
    public:
        double fit(double probabilities[], int num_of_classes, int input_size){
            double result[num_of_classes];
            for (int i = 0; i < num_of_classes; i++){
                double summ = 0;
                for (int i = 0; i < input_size; i++){
                    summ += exp(probabilities[i]);
                };
                result[i] = exp(probabilities[i]) / summ;
            };
            return result[0];
        };
};

class Layer{

    
    public:
        Neuron layer[256];
        int input_size;
        Layer(int input_size){
            for (int i = 0; i < input_size; i++){
                Neuron nn;
                this->layer[i] = nn;
            };
            this->input_size = input_size;
        };
        
};

int main() {
    
    double x[2][6] = {{0.5, 0.2, 0.3, 0.5, 0.7, 0.5}, {0.01, 0.01, 0.001, 0.001, 0.0001, 0.1}};
    double y[2][1] = {{1}, {0}};
    double arrSize = sizeof(y)/sizeof(y[0]);
    
    int layer_1_size = 6;
    int layer_2_size = 3;
    int layer_3_size = 1;

    Layer layer_1(layer_1_size);
    Layer layer_2(layer_2_size);
    Layer layer_3(layer_3_size);

    for (int i = 0; i < layer_1_size; i++){
            double weights[layer_1_size];
            double intercept = (double) rand() / (RAND_MAX);
            for (int k = 0; k < layer_1_size; k++){
                weights[k] = (double) rand() / (RAND_MAX);
                
            };
            layer_1.layer[i].setParams(weights, intercept, layer_1_size);
    };
    
    for (int i = 0; i < layer_2_size; i++){
            double weights[layer_1_size];
            double intercept = (double) rand() / (RAND_MAX);
            for (int k = 0; k < layer_2_size; k++){
                weights[k] = (double) rand() / (RAND_MAX);
                
            };
            layer_2.layer[i].setParams(weights, intercept, layer_2_size);
            
    };

    for (int i = 0; i < layer_3_size; i++){
            double weights[layer_2_size];
            double intercept = (double) rand() / (RAND_MAX);
            for (int k = 0; k < layer_2_size; k++){
                weights[k] = (double) rand() / (RAND_MAX);
                
            }
            layer_3.layer[i].setParams(weights, intercept, layer_3_size);
            
        };

    for (int epoch = 0; epoch < 5; epoch++){
        double lr = 0.1;
        
        for (int i = 0; i < layer_1_size; i++){
            
            layer_1.layer[i].fit(x);
        };

        double output_layer_1[layer_1_size];
        for (int i = 0; i < layer_1_size; i++){
            Sigmoid sigm;
            output_layer_1[i] = sigm.fit(layer_1.layer[i].last_output);
        };


        for (int i = 0; i < layer_2_size; i++){
            
            layer_2.layer[i].fit(output_layer_1);
        };

        double output_layer_2[layer_2_size];
        for (int i = 0; i < layer_2_size; i++){
            Sigmoid sigm;
            output_layer_2[i] = sigm.fit(layer_2.layer[i].last_output);
        };

        for (int i = 0; i < layer_3_size; i++){
            
            layer_2.layer[i].fit(output_layer_2);
        };


        double final_output[layer_3_size];
        for (int i = 0; i < layer_3_size; i++){
            Sigmoid sigm;
            final_output[i] = sigm.fit(layer_3.layer[i].last_output);
            
        };
        cout <<  layer_3.layer[0].weights[0] << endl;
        CrossEntropy cross_before;
        cross_before.setParams(final_output, y, 1);

        for (int k = 0; k < layer_3_size; k++){
            for (int i = 0; i < layer_2_size; i++){
                double x_3 = output_layer_2[i];
                double weight_3 = layer_3.layer[k].weights[i];
                double sigmoid_output_3 = final_output[k];
                double layer_output_3 = layer_3.layer[k].last_output;
                double real = y[0];

                double penality_3 = ((2*real - 1) / (log(2) * sigmoid_output_3)) * 
                (exp(layer_output_3) / (pow((exp(layer_output_3) + 1), 2))) * x_3; 
                
                layer_3.layer[k].weights[i] -= lr * penality_3;

                for (int b = 0; b < layer_2_size;b++){
                    for (int c = 0; c < layer_1_size; c++){
                        double x_2 = output_layer_1[c];
                        double weight_2 = layer_2.layer[b].weights[c];
                        double sigmoid_output_2 = output_layer_2[c];
                        double layer_output_2 = layer_2.layer[b].last_output;
                        double penality_2 = ((2*real - 1) / (log(2) * sigmoid_output_3)) * 
                        (exp(layer_output_3) / (pow((exp(layer_output_3) + 1), 2))) * weight_3 * 
                        (exp(layer_output_2) / (pow((exp(layer_output_2) + 1), 2))) * x_2; 
                        
                        layer_2.layer[b].weights[c] -= lr * penality_2;
                        for (int u = 0; u < layer_1_size; u++){
                            for (int t = 0; t < layer_1_size; t++){
                                double x_1 = x[t];
                                double sigmoid_output_1 = output_layer_1[t];
                                double layer_output_1 = layer_1.layer[u].last_output;
                                double penality_1 = ((2*real - 1) / (log(2) * sigmoid_output_3)) * 
                                (exp(layer_output_3) / (pow((exp(layer_output_3) + 1), 2))) * weight_3 * 
                                (exp(layer_output_2) / (pow((exp(layer_output_2) + 1), 2))) * weight_2 *
                                (exp(layer_output_1) / (pow((exp(layer_output_1) + 1), 2))) * x_1; 
                                
                                layer_1.layer[b].weights[c] -= lr * penality_1;
                            };
                        };
                    };
                };
            }; 
        
        };

       
            
               

               
    }; 
        
};
        
        
        
