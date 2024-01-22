# MagPatRec
Patter Recognition App - Myokinetic Interface

--> Creating data structure (tab)
	* Enter the number of repetitions
	* Choose the movements that you will be loading (order in which you select the movements has to be the same in which you will select the files)
	* Additional: you can add the date of the data in case you want to save a data structure for future use
	* Press on "Load .txt files from a folder" button and choose the txt files of the movements that you selected (be careful to selected in the same order as you selected them in the list box)
	* If there is a separate file for Rest class, choose the file by pressing on the "Load rest data" button

At any point you can save a data structure by going to Data in the menu and choosing "Save magData" 

--> Extracting Features (tab)
	* if you loaded txt files, you can press on "Display magData info" and check if all the information is correct
	* if you didn't load txt files, but you want to load directly the data structure, press on the "Load magData" button
	* In "Magnet features" panel
		* Add "Rest" as movement check box will extract the rest phase between contractions of all movements 
		* You can choose which displacement you want (linear and/or angular) by selecting the checkboxes
		* In case you want to use displacement from the rest position, instead of magnet pair distance, make sure that "Use displacement" check box is on
		* In "Features" list box choose which features you want to extract
		* "Use just samples" check box will not extract features or window them, but use raw data for classification
		* "Sample frequency" will be calculated automatically if it is set as 0
		* "Effective contraction" is a number between 0 and 1, which will define how much of contraction and rest part will be taken into consideration
		* After selecting everything you want, press the "Extract Features"

--> Data Visualization (tab)
	* Press the "Check Data Balance" to see how much data you have for each class
	* If you want to balance the data press on the "Balance data" button
	* You can plot the linear or angular distances, by checking the check boxes and pressing "Plot" button
	* "Check variance" button will give variance of each magnet pair
	* "Plot 3D" is for visualizing data in 3D, select the Algorithm and the feature you want to use for feature reduction

--> Training and Testing (tab)
	* In "Available movements and extracted features" panel choose everything that you want to be considered during classification
	* In "Classifier parameters" panel, choose a classifier you want, the kernel, normalization, and topology
		* If in the topology it is written "test full", it will test on all test data (including transient)
			IMPORTANT: make sure that tmabsTr feature is extracted, it doesn't have to be selected for classification, but has to be extracted
		* "Randomize data" check box will shuffle data within its class (not recommended for cross-validation)
	* Feature selection:
		* The "Percentage of score" box how much of the information do you want to keep
		* In case of Correlation feature selection the "Percentage of score" box is not used, a separate pop-up box will be shown after you press the "Run Offline Training"			
	* To train and test the classifier press the "Run Offline Training" button

--> Training the model for embedded implementation
	* In training and testing tab, if you select topology "Train model", it will train the machine learning model on the entire data set
	* It will additionally show you results of training on a specified percentage of data and testing on rest of data, however the save model is trained on the entire data set
	* In the menu, if you select "Save patRec.txt", it will prepare a .txt file with all necessary parameters for loading into an embedded system

--> Testing model on new data
	* It is possible to load a model and test on new data
	* After you extract the features of your desired test data (as explianed above), you can load a already trained model by clicking in the menu bar "Load patRec"
	* After the model is loaded, choose the "Test model" topology and press the "Run Offline Training" button
