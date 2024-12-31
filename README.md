
# RPI-MD: Protein-RNA interaction prediction by integrating multimodal data via deep learning approaches

In this study, we employed GCN, CNN, FCN, and Attention mechanism to predict interactions between RNA and proteins by collecting severl feature informations, such as sequennce, structure, motif, and physio-chemical.

## Data Organization

The dataset is organized in the following structure:

- `data/`: Contains two subfolders:
  - `raw_data/`: Stores the original dataset.
  - `generated_data/`: Contains data generated for our method.
The program is written in **Python 3.7**. To execute the code, first install the required dependencies by running the following command in your command line:

pip install -r requirements.txt and run the RPI-MD.py method

### Available Methods and Parameters

The following table outlines the available methods and their corresponding parameters:

| Parameter          | Optional Values                               | Description                                                         |
|-------------------|-----------------------------------------------|---------------------------------------------------------------------|
| **Method Name**   | single_dataset_prediction                     | Predict RNA-protein interactions on a chosen dataset.                 |
|                   | compare_different_layers                      | Compare on different number of CNN and GCN layer performance on multiple datasets.                   |
|                   | compare_negative_sample_methods               | Compare different negative sample generation methods on multiple datasets.|
|                   | timeAnalysis                                  | Calculate running time on different datasets.                         |
| **Dataset Name**  | RPI369, RPI2241, RPI1807, NPInter_10412, RPI2241 | Dataset to be used for analysis.                                      |
| **Negative Generation Name** | random, sort_random, sort | Refer to our paper for the specific meaning of these parameters. |
| **Number of Layers** | 1, 2, 3, 4                                   | Number of GCN layers.     
| **Number of Layers** | 1, 2, 3, 4                                   | Number of CNN layers.                                         |

## Hyper-parameters of RPI-MD
We used different hyper-parameters on different datasets and listed them as follows:
|Dataset|Dropout ratio  |Initial learning rate|Weight decay|step_size|γ|Epochs|
|--|--|--|--|--|--|--|
|RPI2241|	0.1|	0.0007|	0.001|	50|	0.7|	25|
|RPI2241_random|	0.1|	0.0007|	0.001|	50|	0.7|	25|
|RPI2241_sort_random|	0.1|	0.0007|	0.001|	50|	0.7|	25|
|RPI369|	0.3|	0.0007|	0.001|	50|	0.7|	50|
|RPI369_random|	0.8|	0.0007|	0.07|	50|	0.7|	50|
RPI369_sort_random|	0.3|	0.0007|	0.001|	50|	0.7|	50|
NPInter10412|	0.1|	0.001|	0.001|	10|	0.7|	30|
NPInter10412_random|	0.1|	0.001|	0.001|	10|	0.7|	30|
NPInter10412_sort_random|	0.1|	0.001|	0.001|	10|	0.7|	30|
RPI2241|	0.1|	0.003|	0.001|	30|	0.7|	50|
RPI2241_random|	0.1|	0.003|	0.001|	30|	0.7|	50|
RPI2241_sort_random|	0.1|	0.003|	0.001|	30|	0.7|	50|
|RPI488|	0.1|	0.0007|	0.001|	50|	0.7|	25|
|RPI488_random|	0.1|	0.0007|	0.001|	50|	0.7|	25|
|RPI488_sort_random|	0.1|	0.0007|	0.001|	50|	0.7|	25|
