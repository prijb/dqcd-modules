#include "DQCD/Modules/interface/BDTinference.h"

BDTinference::BDTinference (std::string filename) {
    XGBoosterCreate(NULL, 0, &booster_);
    XGBoosterLoadModel(booster_, filename.c_str()); // second argument should be a const char *.

    //bst_ulong out_len;
    //const char **out_dump_array;
    //XGBoosterDumpModel(booster_, "", 0, &out_len, &out_dump_array);
    //std::cout << "Dumping" << std::endl;
    //FILE *file = fopen("filename", "w");
    //for (int i=0; i < (int) out_len; i++)
        //fputs(out_dump_array[i], file);
    //fclose(file);
}

// Destructor
BDTinference::~BDTinference() {}

std::vector<float> BDTinference::get_bdt_outputs(std::vector<float> inputs) {
    //bst_ulong num_of_features = 0;
    //XGBoosterGetNumFeature(booster_, &num_of_features);
    //std::cout << "Number of required features: " << num_of_features << std::endl;

    float values[1][inputs.size()];
    int ivar=0;

    //char const config[] =
        //"{\"training\": false, \"type\": 0, "
        //"\"iteration_begin\": 0, \"iteration_end\": 0, \"strict_shape\": false}";
    
    //std::cout << "Size: " << inputs.size() << std::endl;
    for(auto& var : inputs)
    {
        //std::cout << var << std::endl;
        values[0][ivar] = var;
        ++ivar;
    }
    DMatrixHandle dvalues;
    XGDMatrixCreateFromMat(reinterpret_cast<float*>(values), 1, inputs.size(), 0., &dvalues);

    // Dimension of output prediction
    bst_ulong out_dim;

    float const* out_result = NULL;

    auto ret = XGBoosterPredict(
        booster_, dvalues, 0, 0, 0, &out_dim, &out_result);

    XGDMatrixFree(dvalues);

    std::vector<float> results;
    if(ret==0)
    {
        for(unsigned int ic=0; ic < out_dim; ++ic) {
            //std::cout << "Score: " << out_result[ic] << std::endl;
            results.push_back(out_result[ic]);
        }
    }
    //for (auto &elem: results)
        //std::cout << elem << std::endl;

    return results;
}
