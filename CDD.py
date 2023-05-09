from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
from padelpy import from_smiles, padeldescriptor
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm
import logging



#Convert IC50 into pIC50
def pIC50(input):
    pIC50 = []

    for i in input['IC50[nM] Normalized']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    
        
    return input

#Normalaize value of IC50
def norm_value(input):
    norm = []

    for i in input['IC50[nM]']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['IC50[nM] Normalized'] = norm
    
        
    return input

#Padel settings
padeldescriptor(detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=2,
                removesalt=True,
                fingerprints=True)


#Main Class of COMPUTATIONAL DRUG DISCOVERY
class CDD():

    def __init__(self):
        self.data = None
        self.target = new_client.target
        self.activity = new_client.activity
        self.target_selected = None
        self.model = RandomForestRegressor(n_estimators=200)
        self.selector = VarianceThreshold()
        
        
    #Show list of targets from CHEMBL Database
    def show_targets(self, gene: str):
        list_of_targets = pd.DataFrame(self.target.search(gene)).iloc[:, [1,2,5]]
        return list_of_targets
    
    #Download list of molecues from CHEMBL Database for selected target, clean data and load to model
    def select_target(self, gene: str=None, index: int=None, chembl_id: str=None):
        if chembl_id == None:
          self.target_selected = pd.DataFrame(self.target.search(gene))['target_chembl_id'][index]
        else:
          self.target_selected = chembl_id
        activities = self.activity.filter(target_chembl_id=self.target_selected).filter(standard_type="IC50")
        activities = pd.DataFrame.from_dict(activities)
        self.data = activities[['molecule_chembl_id','canonical_smiles','standard_value']]
        self.data = self.data.rename(columns={'standard_value': 'IC50[nM]'})
        #Data cleaning
        self.data = self.data[self.data['IC50[nM]'].notna()]
        self.data = self.data[self.data['canonical_smiles'].notna()]
        self.data = self.data.drop_duplicates(subset=['canonical_smiles'])
        self.data['IC50[nM]'] = self.data['IC50[nM]'].astype('float')
        bioactivity_class = []
        for x in self.data['IC50[nM]']:  
           if x < 1000:
              bioactivity_class.append("Active")
           else:
              bioactivity_class.append("Inactive")
        self.data["bioactivity_class"] = bioactivity_class
        self.data = self.data.reset_index()
        self.data = norm_value(self.data)
        self.data = pIC50(self.data)

    
    #Describe list of molecules smiles using fingerprints from PadelPy. Save result in .CSV file.
    def describe(self):
        smiles = [x for x in self.data['canonical_smiles']]
        # Divide smiles in 10n parts
        if len(smiles) <= 10:
            fingerprints = pd.DataFrame(from_smiles(smiles, fingerprints=True, descriptors=False))
        else:
            n = len(smiles) // 10
            fingerprints_list = []
            for i in tqdm(range(n), desc='Generating fingerprints', unit='part'):
                temp_fingerprints = pd.DataFrame(from_smiles(smiles[i*100:(i+1)*100], fingerprints=True, descriptors=False))
                fingerprints_list.append(temp_fingerprints)
            temp_fingerprints = pd.DataFrame(from_smiles(smiles[n*100:], fingerprints=True, descriptors=False))
            fingerprints_list.append(temp_fingerprints)
            fingerprints = pd.concat(fingerprints_list, ignore_index=True)

        self.data_fingerprints = pd.concat([self.data, fingerprints], axis=1)
        self.data_fingerprints.to_csv(f'{self.target_selected}_fingerprints.csv')
    
    #Train Randomforest model using fingerprints of molecules and pIC50
    #Return list [Score of model, Size of databasa in number of molecues]
    def fit(self):
        try:
            self.data = pd.read_csv(f'{self.target_selected}_fingerprints.csv')
        except FileNotFoundError:
            logging.error('File not found. Please ensure the file exists or use method .describe()')
            return
        x = self.data.iloc[:, 8:]
        y = self.data['pIC50']
        x = self.selector.fit_transform(x)
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.1)
        self.model.fit(x_train, y_train)
        score = self.model.score(x_test, y_test)
        size = len(x_train)
        return [score, size]

    #Predict pIC50 of molecue using model. Return predicted pIC50   
    def predict(self, chembl_id: str):
        #Extract smile of molecue from CHEMBL Database
        molecule = new_client.molecule
        m = pd.DataFrame(molecule.filter(chembl_id=chembl_id))
        m = dict(m.iloc[0]['molecule_structures'])
        smile = [m['canonical_smiles']]
        X = pd.DataFrame(from_smiles(smile, fingerprints=True, descriptors=False))
        X = self.selector.transform(X)
        prediction = self.model.predict(X)
        return prediction[0]

        
        


           
 