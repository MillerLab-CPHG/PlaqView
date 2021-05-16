# PlaqView Data Processing Steps

This is a step-by-step guide on processing scRNA-seq data from public respository. Also serves as a guide for future PlaqView maintainers and for peer-review. This is very detailed and will be continuously updated as needed.

Wei Feng Ma, MSTP UVA. Miller Lab.

## Gathering Data
Currently, we do not have a standard way to deposit and retrieve datasets. Therefore, custom reformatting of each dataset that comes through is required.

For example, Wirka et al., deposited their data in GEO as .txt matrix files that contains all the samples mixed together, whereas Pan et al., had their dataset separated by sample and had to be integrated.

To organize each dataset, the master excel file is edited in the PlaqView directory('Available_dataset.xlsx'). Then, the custom 'DataID' field is used to create a subfolder in the 'data' directory.

The raw input data is not stored within PlaqView project folder, but separately. The preocessed .rds is deposited in the 'data' subdirectory to allow PlaqView to retrieve it.

## Inputing Data
From GEO/Text: copy and run the 'PlaqView_to_seurat_from_GEO.R'.

From Seurat:
Copy and run "PlaqView_preprocessing_master_HUMAN_05-04-2021.R"
Note that this is continuously being updated so please refer to the closest available .r processing file on GitHub!
