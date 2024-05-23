To run :
1. Open IgorPro.
2. Use the file Menu to open procedure file: `HortonBurroneCellReports_0624.ipf`
3. In the Macros Menu, you should now have 2 Sub-Menus: `Load Data` and `Analysis`

4. **To Load the Example Data:**
5. Select the `Analyse all files in folder` Option in the `Load Data` Sub-Menu
   - It will ask you to set the Z step for the confocal sections. For this dataset, leave it at 0.5 microns.
   - Press `Continue` button.
   - Select Folder Containing Data. It is Called `ExampleDataset` 
5. You should now have 13 folders Containing the Data in IgorPro's Data Browser.

6. **To normalize fluorescence per cell:**
   - Select the Folders containing data from the cell to be normalized. See here for instructions on how to work  with the data browser see: [IgorPro Data Browser Tutorial](https://www.wavemetrics.com/news/get-know-feature-data-browser)
   - Select the `Normalize fluor per cell` Option in the `Load Data` Sub-Menu

7. **To Display dignostic results for a single branch** (either inhibitory or excitatory synapses separetely)
   - Place the Red Arrow next to the folder you want to check
   - Select the `Check Results of Branch Analysis` Option in the `Load Data` Sub-Menu

7. **To Analyze the whole dataset at the whole branch level** (As in Fig 1)
   - Select all the Folders you want to include in the analysis in the Data Browser
   - Select the `Summarize Branch Parameters` Option in the `Analysis` Sub-Menu

8. A new folder named `Compile` will be created.
9. It will contain tables called:
    `Shaft_BrSumm` for Inhibitory synapse results and `Spine_BrSumm` for Excitatory synapse results
See here for instructions on how to work with Tables in IgorPro. [Tables in IgorPro](https://www.wavemetrics.com/products/igorpro/creatinggraphs/tables)

These procedures were created and tested in IgorPro versions 9.05 and 6.37.
 


