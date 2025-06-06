# import pacakages
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import requests
import pygmm
from pygmm.model import Scenario
import pygmm.boore_stewart_seyhan_atkinson_2014 as bssa
from scipy.interpolate import interp1d
from scipy.stats import gmean
import statsmodels.api as sm
from scipy.stats import norm
import os

class SSFAdjustment:
    # Constructor
    def __init__(self,T,Longitude,Latitude,Vs30,Mechanism,Region,Record_Set,a,b,Sac_Form):
        self.__T=T # fundamental period of the building
        self.__Longitude=Longitude # Longitude of the site
        self.__Latitude=Latitude # Latitude of the site
        self.__Vs30=Vs30 # Vs30 of the site
        self.__Mechanism=Mechanism # Mechanism of the site
        self.__Region=Region # Region of the site
        self.__Record_Set=Record_Set # NF FF or combined
        self.__Return_Period=9950 # Return Period (9950 should be constant for CMS calculation)
        self.__Z10=-1 # not specified
        self.__Tlist=np.arange(0.01,10.01,0.01) # For calculating CMS
        self.__a=a # lower bound for calculating SaRatio
        self.__b=b # upper bound for calculating SaRatio
        self.__Sac_Form=Sac_Form # Sa collapse Form
        
        # get the Sa_Form and EQ_info Form
        self.__Sa_Form=pd.read_excel(os.path.join(os.path.dirname(__file__),"input_Combined_unscaled.xlsx"))
        self.__EQ_info=pd.read_excel(os.path.join(os.path.dirname(__file__),"EQ_info.xlsx"))

        # Get the corresponding Sa_Form
        if self.__Record_Set == 'Far Field':
            self.__records_spectra = self.__Sa_Form.iloc[:,:45]
        if self.__Record_Set == 'Near Field':
            self.__records_spectra = self.__Sa_Form.iloc[:,[0]+list(range(45,101))]
        if self.__Record_Set == 'Combined':
            self.__records_spectra = self.__Sa_Form

    # Private Method
    def __UHSfromUSGS(self,Longitude,Latitude,Vs30,Return_Period):
        """
        This function can get the UHS from USGS Hazard Toolbox
        inputs:
            Longitude = Longitude of the site [degree] (double)
            Latitude = Latitude of the site [degree] (double)
            Vs30 = Vs30 of the site [m/s] (double)
            Return Period = Return Period of the UHS
        outputs:
        output = a DataFrame includes four columns of data
                1st column: Period [sec]
                2nd column: Sa (the name of the column is the return period) [g]
                3rd column: Magnitude
                4th column: Closest distance [km]
        The function will return 0 if an error happen

        """
        # Check the longitude and latitude is in the range that USGS can accept
        if Longitude <-125 or Longitude>-65:
            raise ValueError("Longitude out of range(bound: [-125,-65])")
            return 0
        if Latitude <24.4 or Latitude>50:
            raise ValueError("LLatitude out of range(bound: [24.4,50])")
            return 0

        # create the url for the request (can cause problem for older version of python, in that case, switch to .format())
        url=f"https://earthquake.usgs.gov/ws/nshmp/conus-2018/dynamic/disagg/{Longitude}/{Latitude}/{Vs30}/{Return_Period}?imt=SA0P01&imt=SA0P02&imt=SA0P03&imt=SA0P05&imt=SA0P075&imt=SA0P1&imt=SA0P15&imt=SA0P2&imt=SA0P25&imt=SA0P3&imt=SA0P4&imt=SA0P5&imt=SA0P75&imt=SA1P0&imt=SA1P5&imt=SA2P0&imt=SA3P0&imt=SA4P0&imt=SA5P0&imt=SA7P5&imt=SA10P0"

        # request for the data
        response=requests.get(url)

        # check the response status
        if response.status_code == 200:
            data=response.json()
        else:
            raise RuntimeError("Fail to get the response from USGS Hazard Toolbox")

        # Build the period list (from the USGS Hazard Tool)
        Period=[0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0]

        # Get the Sa, M and R from the response
        Sa_List=[]
        M_List=[]
        R_List=[]
        for i in range(len(Period)):
            Sa_List.append(data['response']["disaggs"][i]["data"][0]["summary"][0]['data'][2]["value"])
            M_List.append(data['response']["disaggs"][i]["data"][0]['summary'][3]['data'][0]['value'])
            R_List.append(data['response']["disaggs"][i]["data"][0]['summary'][3]['data'][1]['value'])

        # Put all the reulst in a DataFrame
        output_dict = {"T":Period,"Sa (RP={})".format(Return_Period):Sa_List,"M":M_List,"R":R_List}
        output=pd.DataFrame(output_dict)
        return output
    
    def __BSSA2014(self,T, M, Rjb, Vs30, mechanism, region, Z10):
        """
        Calculate the spectral acceleration (median and natural log of the standard deviation)

        Parameters:
        T (float): Period for which to calculate spectral acceleration (in seconds)
        M (float): Moment magnitude of the earthquake (Mw)
        Rjb (float): Joyner-Boore distance to the rupture plane (km)
        Vs30 (float): Time-averaged shear-wave velocity over the top 30m of the site (m/s)
        Z10 (float): Depth to 1.0km/s shear-wave velocity horizon [optional, set to -1 if unspecified] (km)
        mechanism (string): Fault mechanism. Valid options: "U" (unspecified), "SS" (strike-slip),
                            "NS" (normal slip), "RS" (reverse slip)
        region (string): Region for distance attenuation and basin models. Valid options:
                            “global”, “california”, “china”, “italy”, “japan”, “new_zealand”,
                            “taiwan”, “turkey”

        Returns:
        tuple: A tuple containing:
            - numpy.ndarray: The median spectral acceleration for each period
            - numpy.ndarray: The natural log of the standard deviation of spectral acceleration for each period

        """

        # Construct the earthquake scenario
        if Z10 != -1:
            earthquake_scenario = Scenario(mag=M, dist_jb=Rjb, v_s30=Vs30, mechanism=mechanism, region=region, depth_1_0=Z10)
        else:
            earthquake_scenario = Scenario(mag=M, dist_jb=Rjb, v_s30=Vs30, mechanism=mechanism, region=region)

        # Initialize the BSSA2014 model with the scenario
        gmm = bssa.BooreStewartSeyhanAtkinson2014(scenario = earthquake_scenario)

        # Get the median spectral acceleration at the specified period
        median_Sa = gmm.interp_spec_accels(T)

        # Get the natural logarithm of the standard deviation of spectral acceleration
        ln_stdev_Sa = gmm.interp_ln_stds(T, kind='linear')

        return median_Sa, ln_stdev_Sa
    
    def __ComputeTargetEpsilon(self,T, UHS_df, Vs30, region, mechanism, Z10):
        """
        Compute the target epsilon by comparing the GMM spectra to the UHS spectra

        Parameters:
        T (float): Natural period of the building (sec)
        UHS_df (dataframe): Dataframe of UHS spectra and corresponding magnitude and distance for specific site
        Vs30 (float): Time-averaged shear-wave velocity over the top 30m of the site (m/s)
        Z10 (float): Depth to 1.0km/s shear-wave velocity horizon [optional, set to -1 if unspecified] (km)
        mechanism (string): Fault mechanism. Valid options: "U" (unspecified), "SS" (strike-slip),
                            "NS" (normal slip), "RS" (reverse slip)
        region (string): Region for distance attenuation and basin models. Valid options:
                            “global”, “california”, “china”, “italy”, “japan”, “new_zealand”,
                            “taiwan”, “turkey”

        Returns:
        target_epsilon (float): The number of standard deviations the GMM spectra is above the UHS spectra
                                at the building period T
        """

        # Interpolate Sa from the UHS at the building period
        Sa_UHS = np.interp(T, UHS_df['T'], UHS_df.iloc[:,1])

        # Interpolate GMM inputs from the UHS at the building period
        M = np.interp(T, UHS_df['T'], UHS_df['M'])
        Rjb = np.interp(T, UHS_df['T'], UHS_df['R'])

        # Call GMM to calculate median and stdev of Sa(T)
        median_Sa_gmm, ln_stdev_Sa_gmm = self.__BSSA2014(T, M, Rjb, Vs30, mechanism, region, Z10)

        # Interpolate after running GMM - try a few different ways
        # median_Sa_gmm = np.zeros(len(UHS_df))
        # ln_stdev_Sa_gmm = np.zeros(len(UHS_df))

        # for i in len(UHS_df):
        #     T = UHS_df['T'][i]
        #     M = UHS_df['M'][i]
        #     Rjb = UHS_df['R'][i]
        #     median_Sa_gmm, ln_stdev_Sa_gmm = BSSA2014(T, M, Rjb, Vs30, mechanism, region)

        # Compute target epsilon
        target_epsilon = (np.log(Sa_UHS) - np.log(median_Sa_gmm)) / ln_stdev_Sa_gmm

        return target_epsilon
    
    def __compute_rho_BakerJayaram2008(self,Tlist,T):
        """
        Compute the correlation (rho) between spectral accelerations at multiple periods per Baker & Jayaram (2008)

        Parameters:
        Tlist (numpy array): Vector of periods over which the conditional mean spectra will be calculated. Generally, 0.01:0.01:10
        T (float): Natural period of the building (sec)

        Returns:
        rho (numpy array): Vector of correlation factors corresponding to each period in Tlist
        """

        # Initialize output vector
        rho = np.zeros(len(Tlist))

        for idx, period in enumerate(Tlist):

            # Determine Tmin and Tmax
            Tmin = min(period, T)
            Tmax = max(period, T)

            # Calculate coefficient C1
            C1 = 1 - np.cos(np.pi / 2 - 0.366 * np.log(Tmax / max(Tmin, 0.109)))

            # Calculate coefficient C2
            if Tmax < 0.2:
                C2 = 1 - 0.105 * (1 - 1 / (1 + np.exp(100 * Tmax - 5))) * ((Tmax - Tmin) / (Tmax - 0.0099))
            else:
                C2 = 0

            # Calculate coefficient C3
            if Tmax < 0.109:
                C3 = C2
            else:
                C3 = C1

            # Calculate coefficient C4
            C4 = C1 + 0.5 * (np.sqrt(C3) - C3) * (1 + np.cos(np.pi * Tmin / 0.109))

            # Calculate rho
            if Tmax < 0.109:
                rho[idx] = C2
            elif Tmin > 0.109:
                rho[idx] = C1
            elif Tmax < 0.2:
                rho[idx] = min(C2, C4)
            else:
                rho[idx] = C4

        return rho
    
    def __compute_conditional_mean_spectrum(self,Tlist, T, UHS_df, Vs30, mechanism, region, Z10, target_epsilon):
        """
        Compute the conditional mean spectrum (CMS) of a specific building in a specific site

        Parameters:
        Tlist (numpy array): Vector of periods over which the CMS will be calculated. Generally, 0.01:0.01:10
        T (float): Natural period of the building (sec)
        UHS_df (dataframe): Dataframe of UHS spectra and corresponding magnitude and distance for specific site
        Vs30 (float): Time-averaged shear-wave velocity over the top 30m of the site (m/s)
        Z10 (float): Depth to 1.0km/s shear-wave velocity horizon [optional, set to -1 if unspecified] (km)
        mechanism (string): Fault mechanism. Valid options: "U" (unspecified), "SS" (strike-slip),
                            "NS" (normal slip), "RS" (reverse slip)
        region (string): Region for distance attenuation and basin models. Valid options:
                            “global”, “california”, “china”, “italy”, “japan”, “new_zealand”,
                            “taiwan”, “turkey”
        target_epsilon (float): The number of standard deviations the GMM spectra is above the UHS spectra
                                at the building period T

        Returns:
        median_Sa_CMS (numpy array): The median CMS spectral acceleration at each period in Tlist
        ln_stdev_Sa_CMS (numpy array): The natural log of the standard deviation of the CMS spectral acceleration
                                        at each period in Tlist
        """

        # Compute rho
        rho = self.__compute_rho_BakerJayaram2008(Tlist,T)

        # Interpolate GMM inputs from the UHS at the building period
        M = np.interp(T, UHS_df['T'], UHS_df['M'])
        Rjb = np.interp(T, UHS_df['T'], UHS_df['R'])

        # Compute Sa spectrum using the GMM over the period range Tlist
        median_Sa_gmm = np.zeros(len(Tlist))
        ln_stdev_Sa_gmm = np.zeros(len(Tlist))
        for idx, period in enumerate(Tlist):
            median_Sa_gmm[idx], ln_stdev_Sa_gmm[idx] = self.__BSSA2014(period, M, Rjb, Vs30, mechanism, region, Z10)

        # "Condition" the spectrum with rho and target epsilon
        median_Sa_CMS = np.exp(np.log(median_Sa_gmm) + ln_stdev_Sa_gmm * rho * target_epsilon)
        ln_stdev_Sa_CMS = ln_stdev_Sa_gmm * np.sqrt(1 - rho**2)

        return median_Sa_CMS, ln_stdev_Sa_CMS
    
    def __calculate_target_SaRatio(self, median_Sa_CMS, a, b, T, Tlist):

        # Slice conditional mean spectra from Ta to Tb
        mask = (Tlist >= T*a) & (Tlist <= T*b)
        periods_range = Tlist[mask]
        Sa_CMS_range = median_Sa_CMS[mask]

        # Add in Ta if its not the first period listed in periods_range
        if periods_range[0] != T*a:
            Sa_CMS_Ta = np.interp(T*a, Tlist, median_Sa_CMS)
            Sa_CMS_Ta = np.array([Sa_CMS_Ta])
            Sa_CMS_range = np.concatenate((Sa_CMS_Ta, Sa_CMS_range))

        # Add in Tb if its not the last period listed in periods_range
        if periods_range[-1] != T*b:
            Sa_CMS_Tb = np.interp(T*b, Tlist, median_Sa_CMS)
            Sa_CMS_Tb = np.array([Sa_CMS_Tb])
            Sa_CMS_range = np.concatenate((Sa_CMS_range, Sa_CMS_Tb))

        # Compute Sa at building period T using conditional mean spectra
        Sa_CMS_T = np.interp(T, Tlist, median_Sa_CMS)

        # Compute target SaRatio
        target_SaRatio = Sa_CMS_T / gmean(Sa_CMS_range)

        return target_SaRatio
    
    def __Calculate_epsilon_of_RecordSet(self, T, RecordSet, EQ_info, records_spectra,Z10):

        # Import 5% damped response spectra for each record in the specified RecordSet
        if RecordSet == 'Far Field':
            records_list = EQ_info['EQ_ID'][:44].to_list()
        if RecordSet == 'Near Field':
            records_list = EQ_info['EQ_ID'][44:].to_list()
        if RecordSet == 'Combined':
            records_list = EQ_info['EQ_ID'].to_list()

        # Define periods associated with the record spectra
        periods = records_spectra['T']

        # Compute spectral acceleration at T for both target spectra and record spectra
        median_Sa_gmm = np.zeros(len(records_list))
        ln_stdev_Sa_gmm = np.zeros(len(records_list))
        Sa_records = np.zeros(len(records_list))

        for idx,record in enumerate(records_list):

            # Compute Sa_GMM (SaT1 of target spectra) for each record site via GMM
            M = EQ_info.loc[EQ_info['EQ_ID'] == record, 'M'].values[0]
            Rjb = EQ_info.loc[EQ_info['EQ_ID'] == record, 'Rjb'].values[0]
            mechanism = EQ_info.loc[EQ_info['EQ_ID'] == record, 'mechanism'].values[0]
            Vs30 = EQ_info.loc[EQ_info['EQ_ID'] == record, 'Vs30'].values[0]
            region = EQ_info.loc[EQ_info['EQ_ID'] == record, 'region'].values[0]

            median_Sa_gmm[idx], ln_stdev_Sa_gmm[idx] = self.__BSSA2014(T, M, Rjb, Vs30, mechanism, region, Z10)

            # Compute Sa_record (SaT1 of record spectra) for each record in the set
            record_spectra = records_spectra[record]
            Sa_records[idx] = np.interp(T, periods, record_spectra)

        # Compute epsilon for each record in the set
        epsilon = (np.log(Sa_records) - np.log(median_Sa_gmm)) / ln_stdev_Sa_gmm

        # Prepare output dataframe with EQ_ID and epsilon
        epsilon={"EQ_ID":records_list,"epsilon":epsilon}
        epsilon=pd.DataFrame(epsilon)

        return epsilon
    
    def __SaRatioForRecord(self, Sa_Form,T,Ta_T,Tb_T):
        """
        # This function can get SaRtio for the records
        inputs:
        # Sa_Form = Sa for all the ground motion records (DataFrame)
        # T = Fundamental period of the structure [s] (double)
        # Ta_T = the ratio of Ta to T [no unit] (double)
        # Tb_T = the ratio of Tb to T [no unit] (double)
        outputs:
        # output = a DataFrame includes two columns of data
        #          1st column: EQ_ID (int)
        #          2nd column: Sa Ratio (double)
        # The function will return 0 if an error happen
        """
        # check the fromat of the DataFrame
        Column_Names=Sa_Form.columns
        # check the first column is T
        if Column_Names[0] != "T":
            raise ValueError("The first columne should be the period, and please name the first column 'T'")
            return 0
        # check other columns are Sa for each gm record
        for i in range(1,Sa_Form.shape[1]):
            if type(Column_Names[i]) != int:
                raise ValueError("The column {} should be the Sa for the ground motion record, and please name the column using EQ_ID".format(i+1))
                return 0
        # check the value in the spreadsheet
        Column_Types=Sa_Form.dtypes
        for i in range(len(Column_Types)):
            if Column_Types.iloc[i]!=float:
                raise ValueError("Column {} is not the float number".format(i+1))

        # build the fist col for the output dataframe
        record_name=Sa_Form.columns[1:]
        output_dict={"EQ_ID":record_name}

        # Define Tmin and Tmax
        T_a=Ta_T*T
        T_b=Tb_T*T

        # Loop to calculate the SA Ratio for each record in the set
        Sa_Ratio=[]
        for i in range(1,Sa_Form.shape[1]):
            # Find the SA at T
            for j in range(Sa_Form.shape[0]):
                if Sa_Form.iloc[j,0]>T:
                    Sa_T=np.interp(T,Sa_Form.iloc[j-1:j+1,0],Sa_Form.iloc[j-1:j+1,i])
                    break
            # Find the SA at [T_a:0.01:T_b]
            Sa_average_list=[]
            # Find Sa for Ta in the form
            for j in range(Sa_Form.shape[0]):
                if Sa_Form.iloc[j,0]>T_a:
                    Sa_Ta=np.interp(T_a,Sa_Form.iloc[j-1:j+1,0],Sa_Form.iloc[j-1:j+1,i])
                    Ta_index_start=j
                    break
            # Find Sa for Tb in the form
            for j in range(Sa_Form.shape[0]):
                if Sa_Form.iloc[j,0]>T_b:
                    Sa_Tb=np.interp(T_b,Sa_Form.iloc[j-1:j+1,0],Sa_Form.iloc[j-1:j+1,i])
                    Tb_index_end=j-1
                    break
            Sa_average_list.append(Sa_Ta)
            Sa_average_mid=Sa_Form.iloc[Ta_index_start:Tb_index_end+1,i]
            Sa_average_mid.tolist()
            Sa_average_list.extend(Sa_average_mid)
            Sa_average_list.append(Sa_Tb)
            # Calculate the geomean to get the Sa_average
            Sa_avg=gmean(Sa_average_list)
            # Calculate the SA Ratio
            Sa_Ratio.append(Sa_T/Sa_avg)
        output_dict["Sa Ratio"]=Sa_Ratio
        output=pd.DataFrame(output_dict)
        output=output.sort_values(by="EQ_ID")
        return output
    
    def __RegressionEpsilonSacollapse(self, Sac_Form,Epsilon_Form):
        """
        # This function can regression result of Sac-epsilon relationship in log scale
        inputs:
        # Sac_Form = Sa collapse for all the ground motion records (DataFrame)
        #            1st column: EQ_ID (int)
        #            2nd column: Sa collapse (float)
        # Epsilon_Form = Epsilon at the fundamental period for all the ground motion records
        #            1st column: EQ_ID (int)
        #            2nd column: Epsilon (float)
        outputs:
        # output = a dictionary:
        #          "Original Data": EQ_ID (list),Sac (list),Epsilon (list) from the input DataFrame
        #          "Plot Data": Sac (list) and Epsilon (list) fit for the regression
        #          "beta":beta (list) of the regression
        #          "Correlation":correlation (float) of epsilon and sac
        # The function will return 0 if an error happen
        """
        # rearange the two DataFrame to make sure the EQ_ID are in the same order
        Sac_Form=Sac_Form.sort_values(by=Sac_Form.columns[0])
        Epsilon_Form=Epsilon_Form.sort_values(by=Epsilon_Form.columns[0])
        # check the name of the records corresponds in the two DataFrame
        for i in range(Sac_Form.shape[0]):
            if Sac_Form.iloc[i,0]!=Epsilon_Form.iloc[i,0]:
                raise ValueError("The Ground Motion Records in Sa Collapse and Epsilon are not the same")
                return 0
        # Calculate the relationship between Sa collapse and Epsilon
        X=np.ones([Epsilon_Form.shape[0],2])
        # Loop to fill in the Epsilon
        for i in range(Epsilon_Form.shape[0]):
            X[i][0]=Epsilon_Form.iloc[i,1]
        Y=np.log(Sac_Form.iloc[:,1].tolist())
        beta=np.linalg.inv(X.T@X)@X.T@Y
        x_plot=np.arange(np.min(Epsilon_Form.iloc[:,1]),np.max(Epsilon_Form.iloc[:,1]),0.01)
        y_plot=np.exp(x_plot*beta[0]+beta[1])
        cof=np.corrcoef(Epsilon_Form.iloc[:,1],Y)

        # Build the output dict
        # The Original Data
        Original_Data={"EQ_ID":Sac_Form.iloc[:,0].tolist(),
                    "Sac":Sac_Form.iloc[:,1].tolist(),
                    "Epsilon":Epsilon_Form.iloc[:,1].tolist()}
        # The fitting data using the regression result
        Plot_Data={"Sac":y_plot.tolist(),
                "Epsilon":x_plot.tolist()}
        # Put all the result in one dict
        Output={"Original Data":Original_Data,
            "Plot Data":Plot_Data,
            "beta":beta.tolist(),
            "Correlation":cof[0][1]}

        return Output
    
    def __PlotEpsilonSac(self, Result_File):
        """
        # This function can plot the Epsilon vs Sa Collapse
        inputs:
        # The dictionary from part 8:
        #          "Original Data": EQ_ID (list),Sac (list),Epsilon (list) from the input DataFrame
        #          "Plot Data": Sac (list) and Epsilon (list) fit for the regression
        #          "beta":beta (list) of the regression
        #          "Correlation":correlation (float) of Epsilon and Sac
        outputs:
        # No Output, Only a figure
        """
        # create the plot
        fig,ax=plt.subplots()

        # Plot the original data in scatter points
        ax.scatter(Result_File["Original Data"]["Epsilon"],Result_File["Original Data"]["Sac"])

        # Plot the Regression Curve
        ax.plot(Result_File["Plot Data"]["Epsilon"],Result_File["Plot Data"]["Sac"],label=r"ln(y)={}x+{}, $\rho$={}".format(round(Result_File["beta"][0],2),round(Result_File["beta"][1],2),round(Result_File["Correlation"],2)))

        ax.set_xlabel(r"$\epsilon(T_1)$")
        ax.set_ylabel(r"$Sa_{collapse}$")
        ax.set_title(r"$Sa_{collapse}$-$\epsilon$ Relationship")
        ax.grid(True)
        ax.legend()
        return fig

    def __RegressionSaratioSacollapse(self, Sac_Form,SaRatio_Form):
        """
        # This function can regression result of Sac-Sa Ratio relationship in log scale
        inputs values
        # Sac_Form = Sa collapse for all the ground motion records (DataFrame)
        #            1st column: EQ_ID (int)
        #            2nd column: Sa collapse (float)
        # SaRatio_Form = Sa Ratio at the fundamental period for all the ground motion records
        #            1st column: EQ_ID (int)
        #            2nd column: Sa Ratio (float)
        outputs:
        # output = a dictionary:
        #          "Original Data": EQ_ID (list),Sac (list),SaRatio (list) from the input DataFrame
        #          "Plot Data": Sac (list) and SaRatio (list) fit for the regression
        #          "beta":beta (list) of the regression
        #          "Correlation":correlation (float) of epsilon and sac
        # The function will return 0 if an error happen
        """
        # rearange the two DataFrame to make sure the EQ_ID are in the same order
        Sac_Form=Sac_Form.sort_values(by=Sac_Form.columns[0])
        SaRatio_Form=SaRatio_Form.sort_values(by=SaRatio_Form.columns[0])
        # check the name of the records corresponds in the two DataFrame
        for i in range(Sac_Form.shape[0]):
            if Sac_Form.iloc[i,0]!=SaRatio_Form.iloc[i,0]:
                raise ValueError("The Ground Motion Records in Sa Collapse and Sa Ratio are not the same")
                return 0
        # Calculate the relationship between Sa collapse and Epsilon
        X=np.ones([SaRatio_Form.shape[0],2])
        # Loop to fill in the Epsilon
        for i in range(SaRatio_Form.shape[0]):
            X[i][0]=SaRatio_Form.iloc[i,1]
        Y=Sac_Form.iloc[:,1].tolist()
        beta=np.linalg.inv(X.T@X)@X.T@Y
        x_plot=np.arange(np.min(SaRatio_Form.iloc[:,1]),np.max(SaRatio_Form.iloc[:,1]),0.01)
        y_plot=x_plot*beta[0]+beta[1]
        cof=np.corrcoef(SaRatio_Form.iloc[:,1],Y)

        # Build the output dict
        # The Original Data
        Original_Data={"EQ_ID":Sac_Form.iloc[:,0].tolist(),
                    "Sac":Sac_Form.iloc[:,1].tolist(),
                    "Sa Ratio":SaRatio_Form.iloc[:,1].tolist()}
        # The fitting data using the regression result
        Plot_Data={"Sac":y_plot.tolist(),
                "Sa Ratio":x_plot.tolist()}
        # Put all the result in one dict
        Output={"Original Data":Original_Data,
            "Plot Data":Plot_Data,
            "beta":beta.tolist(),
            "Correlation":cof[0][1]}

        return Output
    
    def __PlotSaratioSac(self, Result_File):
        """
        # This function can plot the Sa Ratio vs Sa Collapse
        inputs:
        # The dictionary from part 10:
        #          "Original Data": EQ_ID (list),Sac (list),Sa Ratio (list) from the input DataFrame
        #          "Plot Data": Sac (list) and Sa Ratio (list) fit for the regression
        #          "beta":beta (list) of the regression
        #          "Correlation":correlation (float) of Sa Ratio and Sac
        outputs:
        # No Output, Only a figure
        """
        # create the plot
        fig,ax=plt.subplots()

        # Plot the original data in scatter points
        ax.scatter(Result_File["Original Data"]["Sa Ratio"],Result_File["Original Data"]["Sac"])

        # Plot the Regression Curve
        ax.plot(Result_File["Plot Data"]["Sa Ratio"],Result_File["Plot Data"]["Sac"],label=r"y={}x+{}, $\rho$={}".format(round(Result_File["beta"][0],2),round(Result_File["beta"][1],2),round(Result_File["Correlation"],2)))

        ax.set_xlabel(r"$Sa\ Ratio$")
        ax.set_ylabel(r"$Sa_{collapse}$")
        ax.set_title(r"$Sa_{collapse}$-$Sa\ Ratio$ Relationship")
        ax.grid(True)
        ax.legend()
        return fig

    def __SSFForEpsilon(self, Epsilon_Target,Result_File):
        """
        # This function can calculate the SSF for Epsilon
        inputs:
        # Epsilon_Target = Target Epsilon from step 3 (float)
        # Result_File = The dictionary from step 8:
        #          "Original Data": EQ_ID (list),Sac (list),Epsilon (list) from the input DataFrame
        #          "Plot Data": Sac (list) and Epsilon (list) fit for the regression
        #          "beta":beta (list) of the regression
        #          "Correlation":correlation (float) of Epsilon and Sac
        outputs:
        # SSF = spectra shape factor (float)
        # The function will return 0 if an error happens
        """

        # calculate the mean value of the Epsilon from the dictionary
        Epsilon_Set=np.mean(Result_File["Original Data"]["Epsilon"])

        # calculate the SSF
        SSF=np.exp(Result_File["beta"][0]*(Epsilon_Target-Epsilon_Set))
        return SSF
    
    def __SSFForSaRatio(self, Sa_Ratio_Target,Result_File):
        """
        # This function can calculate the SSF for Sa Ratio
        inputs:
        # Sa_Ratio_Target = Target Sa Ratio from step 5 (float)
        # Result_File = The dictionary from step 10:
        #          "Original Data": EQ_ID (list),Sac (list),Sa Ratio (list) from the input DataFrame
        #          "Plot Data": Sac (list) and Sa Ratio (list) fit for the regression
        #          "beta":beta (list) of the regression
        #          "Correlation":correlation (float) of Sa Ratio and Sac
        outputs:
        # SSF = spectra shape factor (float)
        # The function will return 0 if an error happens
        """
        # calculate the mean value of the Epsilon from the dictionary
        Sa_Ratio_Set=np.mean(Result_File["Original Data"]["Sa Ratio"])

        # calculate the SSF
        SSF=Result_File["beta"][0]*(Sa_Ratio_Target-Sa_Ratio_Set)
        return SSF
    
    def __mle_fragility(self, IM, num_obs, num_fail):
        IM = np.asarray(IM).reshape(-1)
        num_fail = np.asarray(num_fail).reshape(-1)

        if np.isscalar(num_obs):
            num_obs = np.ones_like(num_fail) * num_obs
        else:
            num_obs = np.asarray(num_obs).reshape(-1)

        num_success = num_obs - num_fail

        Y = np.column_stack((num_fail, num_success))

        X = sm.add_constant(np.log(IM))

        model = sm.GLM(Y, X, family=sm.families.Binomial(link=sm.families.links.Probit()))
        results = model.fit()

        b0, b1 = results.params
        theta = np.exp(-b0 / b1)
        beta = 1 / b1

        return theta, beta
    
    def __MedianDispersionForEpsilonCurve(self, Sac_Form,SSF):
        """
        # This function can calculate the median and dispersion for Epsilon shifted fragility curve
        inputs:
        # Sac_Form = Sa collapse for all the ground motion records (DataFrame)
        #            1st column: EQ_ID (int)
        #            2nd column: Sa collapse (float)
        # SSF = Spectra shape factor from step 12 (float)
        outputs:
        # A dictionary includes:
        #    "Type": the date is using "Epsilon" shift (string:"Epsilon")
        #    "Meidan Sa Collapse": The median Sa Collapse of the fragility curve (float)
        #    "Shifted Sa Collapse": The shifted median Sa Collapse of the fragility curve (float)
        #    "Dispersion Sa Collapse": The dispersion of the fragility curve(float)
        # The function will return 0 if an error happens
        """
        # Get the median and the dispersion for the Sac
        Median_Sac_Original,Dispersion_Sac_Original=self.__mle_fragility(sorted(Sac_Form.iloc[:,1]),len(Sac_Form.iloc[:,1]),np.arange(1,Sac_Form.shape[0]+1))

        # Shift the median according to the SSF
        Median_Sac_Shift=SSF*Median_Sac_Original

        # create the output dict
        Output= {
            "Type": "Epsilon",
            "Median Sa Collapse":Median_Sac_Original,
            "Shifted Sa Collapse":Median_Sac_Shift,
            "Dispersion Sa Collapse": Dispersion_Sac_Original
        }
        return Output
    
    def __MedianDispersionForSaRatioCurve(self, Sac_Form,SSF):
        """
        # This function can calculate the median and dispersion for Sa Ratio shifted fragility curve
        inputs:
        # Sac_Form = Sa collapse for all the ground motion records (DataFrame)
        #            1st column: EQ_ID (int)
        #            2nd column: Sa collapse (float)
        # SSF = Spectra shape factor from step 13 (float)
        outputs:
        # A dictionary includes:
        #    "Type": the date is using "Sa Ratio" shift (string:"Sa Ratio")
        #    "Meidan Sa Collapse": The median Sa Collapse of the fragility curve (float)
        #    "Shifted Sa Collapse": The shifted median Sa Collapse of the fragility curve (float)
        #    "Dispersion Sa Collapse": The dispersion of the fragility curve(float)
        # The function will return 0 if an error happens
        """
        # Get the median and the dispersion for the Sac
        Median_Sac_Original,Dispersion_Sac_Original=self.__mle_fragility(sorted(Sac_Form.iloc[:,1]),len(Sac_Form.iloc[:,1]),np.arange(1,Sac_Form.shape[0]+1))

        # Shift the median according to the SSF
        Median_Sac_Shift=SSF+Median_Sac_Original

        # create the output dict
        Output= {
            "Type": "Sa Ratio",
            "Median Sa Collapse":Median_Sac_Original,
            "Shifted Sa Collapse":Median_Sac_Shift,
            "Dispersion Sa Collapse": Dispersion_Sac_Original
        }
        return Output
    
    def __PlotFragilityCurve(self, Paras,Option="both"):
        """
        # This function can plot the fragility curves for different options
        inputs:
        # Paras = A list of dictionaries, acceptable lists are:
        #         [dictionary from step14]
        #         [dictionary from step15]
        #         [dictionary from step14, dictionary from step15]
        # Option = "Epsilon" (string): Only plot the Epsilon shift
        #          "Sa Ratio" (string): Only plot the Sa Ratio shift
        #          "Both" (string, default) : Plot all the curves
        outputs:
        # No Output, Only a figure
        # The function will return 0 if an error happens
        """
        # check the Paras is a list
        if type(Paras) != list:
            raise ValueError("The first input should be a list")
            return 0

        # the "Epsilon" option
        if Option == "Epsilon":
            # check the paras are correct
            for i in range(len(Paras)):
                if Paras[i]["Type"] == "Epsilon":
                    Paras_dict=Paras[i]
                    break
                if i == len(Paras)-1:
                    raise ValueError("No available parameters for Epsilon shift")
                    return 0
            # Calculate x and y for the plot
            x=np.arange(0.01,8,0.01)
            y_Original=norm.cdf((np.log(x)-np.log(Paras_dict["Median Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            y_Epsilon=norm.cdf((np.log(x)-np.log(Paras_dict["Shifted Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            # Plot the fragility curve
            fig,ax=plt.subplots()
            ax.plot(x,y_Original,label=r"IDA, $\theta$={},$\beta$={}".format(round(Paras_dict["Median Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))
            ax.plot(x,y_Epsilon,label=r"$\epsilon$, $\theta$={},$\beta$={}".format(round(Paras_dict["Shifted Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))
            ax.legend()
            ax.grid(True)
            ax.set_xlabel(r"$Sa$")
            ax.set_ylabel(r"P(c)")
            ax.set_title("Fragility Curves")
            return fig
        # the "Sa Ratio" option
        elif Option == "Sa Ratio":
            # check the paras are correct
            for i in range(len(Paras)):
                if Paras[i]["Type"] == "Sa Ratio":
                    Paras_dict=Paras[i]
                    break
                if i == len(Paras)-1:
                    raise ValueError("No available parameters for Sa Ratio shift")
                    return 0
            # Calculate x and y for the plot
            x=np.arange(0.01,8,0.01)
            y_Original=norm.cdf((np.log(x)-np.log(Paras_dict["Median Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            y_Epsilon=norm.cdf((np.log(x)-np.log(Paras_dict["Shifted Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            # Plot the fragility curve
            fig,ax=plt.subplots()
            ax.plot(x,y_Original,label=r"IDA, $\theta$={},$\beta$={}".format(round(Paras_dict["Median Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))
            ax.plot(x,y_Epsilon,label=r"$Sa\ Ratio$, $\theta$={},$\beta$={}".format(round(Paras_dict["Shifted Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))
            ax.legend()
            ax.grid(True)
            ax.set_xlabel(r"$Sa$")
            ax.set_ylabel(r"P(c)")
            ax.set_title("Fragility Curves")
            return fig
        # the "both" option
        elif Option=="both":
            # first plot the Epsilon
            # check the paras are correct
            for i in range(len(Paras)):
                if Paras[i]["Type"] == "Epsilon":
                    Paras_dict=Paras[i]
                    break
                if i == len(Paras)-1:
                    raise ValueError("No available parameters for Epsilon shift")
                    return 0
            # Calculate x and y for the plot
            x=np.arange(0.01,8,0.01)
            y_Original=norm.cdf((np.log(x)-np.log(Paras_dict["Median Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            y_Epsilon=norm.cdf((np.log(x)-np.log(Paras_dict["Shifted Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            # Plot the fragility curve
            fig,ax=plt.subplots()
            ax.plot(x,y_Original,label=r"IDA, $\theta$={},$\beta$={}".format(round(Paras_dict["Median Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))
            ax.plot(x,y_Epsilon,label=r"$\epsilon$, $\theta$={},$\beta$={}".format(round(Paras_dict["Shifted Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))

            # then plot the Sa Ratio
            # check the paras are correct
            for i in range(len(Paras)):
                if Paras[i]["Type"] == "Sa Ratio":
                    Paras_dict=Paras[i]
                    break
                if i == len(Paras)-1:
                    raise ValueError("No available parameters for Sa Ratio shift")
                    return 0
            # Calculate x and y for the plot
            y_Epsilon=norm.cdf((np.log(x)-np.log(Paras_dict["Shifted Sa Collapse"]))/Paras_dict["Dispersion Sa Collapse"])
            # Plot the fragility curve
            ax.plot(x,y_Epsilon,label=r"$Sa\ Ratio$, $\theta$={},$\beta$={}".format(round(Paras_dict["Shifted Sa Collapse"],2),round(Paras_dict["Dispersion Sa Collapse"],2)))
            ax.legend()
            ax.grid(True)
            ax.set_xlabel(r"$Sa$")
            ax.set_ylabel(r"P(c)")
            ax.set_title("Fragility Curves")
            return fig
        else:
            raise ValueError("No such option, availbale options:'Epsilon','Sa Ratio','both'")
            return 0
            
    def __PlotSacT(self, T,Sac_form,Sa_form):
        """
        Input:
            T = the fundamental period of the structure
            Sac_form = Sa collapse for all the ground motion records [g] (DataFrame)
                1st column: EQ_ID (int)
                2nd column: Sa collapse (float)
            Sa_form = Sa for all the ground motion records [g] (DataFrame)
        Output:
            No output,just a figure
        """

        # loop through all the ground motion records
        fig,ax=plt.subplots()
        column_names=Sa_form.columns
        Sac_form=Sac_form.sort_values(by=Sac_form.columns[1])
        Sac_plot=Sac_form.iloc[[0,1,2,-1,-2,-3],:]
        for i in range(Sac_plot.shape[0]):
            # find the right EQ_ID
            for j in range(1,Sa_form.shape[1]):
                if Sac_plot.iloc[i,0] == column_names[j]:
                    # find the Sa(T)
                    for k in range(Sa_form.shape[0]):
                        if Sa_form.iloc[k,0]>T:
                            Sa_T=round(np.interp(T,Sa_form.iloc[k-1:k+1,0],Sa_form.iloc[k-1:k+1,j]),4)
                            SF=Sac_plot.iloc[i,1]/Sa_T
                            ax.plot(Sa_form.iloc[3:,0],Sa_form.iloc[3:,j]*SF)
                            break
                    break
        ymin, ymax = plt.ylim()
        ax.grid(True)
        ax.set_xlabel("Period [s]")
        ax.set_ylabel("Sa")
        ax.plot([T,T],[ymin,ymax],color="black",ls="--")
        ax.set_ylim([0,4])
        ax.set_xlim([0,4])
        ax.set_title("Sa for different periods at Sa collapse")
        return fig
    
    # Public Method
    # Function to run the analysis
    def Run(self):
        self.__UHS_df=self.__UHSfromUSGS(self.__Longitude,self.__Latitude,self.__Vs30,self.__Return_Period)
        self.__target_epsilon = self.__ComputeTargetEpsilon(self.__T, self.__UHS_df, self.__Vs30, self.__Region, self.__Mechanism, self.__Z10)
        self.__median_Sa_CMS,self.__ln_stdev_Sa_CMS = self.__compute_conditional_mean_spectrum(self.__Tlist, self.__T, self.__UHS_df, self.__Vs30, self.__Mechanism, self.__Region, self.__Z10, self.__target_epsilon)
        self.__target_SaRatio = self.__calculate_target_SaRatio(self.__median_Sa_CMS, self.__a, self.__b, self.__T, self.__Tlist)
        self.__epsilon = self.__Calculate_epsilon_of_RecordSet(self.__T, self.__Record_Set, self.__EQ_info, self.__records_spectra,self.__Z10)
        self.__SaRatio = self.__SaRatioForRecord(self.__records_spectra,self.__T,self.__a,self.__b)
        self.__Regression_epsilon = self.__RegressionEpsilonSacollapse(self.__Sac_Form,self.__epsilon)
        # self.__PlotEpsilonSac(self.__Regression_epsilon)
        self.__Regression_SaRatio=self.__RegressionSaratioSacollapse(self.__Sac_Form,self.__SaRatio)
        # self.__PlotSaratioSac(self.__Regression_SaRatio)
        self.__SSF_epsilon=self.__SSFForEpsilon(self.__target_epsilon,self.__Regression_epsilon)
        self.__SSF_SaRatio=self.__SSFForSaRatio(self.__target_SaRatio,self.__Regression_SaRatio)
        self.__Fragility_Curve_epsilon=self.__MedianDispersionForEpsilonCurve(self.__Sac_Form,self.__SSF_epsilon)
        self.__Fragility_Curve_SaRatio=self.__MedianDispersionForSaRatioCurve(self.__Sac_Form,self.__SSF_SaRatio)
        # self.__PlotFragilityCurve([self.__Fragility_Curve_epsilon,self.__Fragility_Curve_SaRatio],Option="both")
        # self.__PlotSacT(self.__T,self.__Sac_Form,self.__records_spectra)
    # Function to plot Sac vs epsilon
    def SacEpsilonPlot(self):
        fig=self.__PlotEpsilonSac(self.__Regression_epsilon)
        return fig
    # Function to plot Sac vs Sa Ratio
    def SacSaRatioPlot(self):
        fig=self.__PlotSaratioSac(self.__Regression_SaRatio)
        return fig
    # Function to plot the fragility curve
    def FragilityCurve(self,option):
        fig=self.__PlotFragilityCurve([self.__Fragility_Curve_epsilon,self.__Fragility_Curve_SaRatio],option)
        return fig
    
    # Function to plot the Sac
    def SacPlot(self):
        fig=self.__PlotSacT(self.__T,self.__Sac_Form,self.__records_spectra)
        return fig 

    @ property
    def SSF_epsilon(self):return self.__Fragility_Curve_epsilon["Shifted Sa Collapse"]/self.__Fragility_Curve_epsilon["Median Sa Collapse"]
    @ property
    def SSF_Sa_Ratio(self):return self.__Fragility_Curve_SaRatio["Shifted Sa Collapse"]/self.__Fragility_Curve_SaRatio["Median Sa Collapse"]
    @ property
    def UHS(self): return self.__UHS_df
    @ property
    def CMS(self):
        Output_CMS={
            "T":self.__Tlist,
            "Meidan":self.__median_Sa_CMS,
            "Ln_std":self.__ln_stdev_Sa_CMS
        }
        Output_CMS=pd.DataFrame(Output_CMS)
        return Output_CMS
    @ property
    def epsilon_for_record(self): return self.__epsilon
    @ property
    def Sa_Ratio_for_record(self): return self.__SaRatio
    @ property
    def epsilon_Regression_Result(self):
        Output_Sac_epsilon_regression={"Sa_collapse":self.__Regression_epsilon['Plot Data']['Sac'],"epsilon":self.__Regression_epsilon['Plot Data']['Epsilon']}
        Output_Sac_epsilon_regression=pd.DataFrame(Output_Sac_epsilon_regression)
        return Output_Sac_epsilon_regression
    @ property
    def Sa_Ratio_Regression_Result(self):
        Output_Sac_SaRatio_regression={"Sa_collapse":self.__Regression_SaRatio['Plot Data']['Sac'],"Sa_Ratio":self.__Regression_SaRatio['Plot Data']['Sa Ratio']}
        Output_Sac_SaRatio_regression=pd.DataFrame(Output_Sac_SaRatio_regression)
        return Output_Sac_SaRatio_regression
    
    # Gain the input data
    @ property
    def T(self):return self.__T
    @ property
    def Longitude(self):return self.__Longitude
    @ property
    def Latitude(self):return self.__Latitude
    @ property
    def Vs30(self):return self.__Vs30
    @ property
    def Mechanism(self):return self.__Mechanism
    @ property
    def Region(self):return self.__Region
    @ property
    def Record_Set(self):return self.__Record_Set
    @ property
    def a(self):return self.__a
    @ property
    def b(self):return self.__b
    @ property
    def SacForm(self):return self.__Sac_Form
