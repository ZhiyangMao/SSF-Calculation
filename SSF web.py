import streamlit as st
from OOP_SSF_Adjustment import SSFAdjustment
import pandas as pd

st.set_page_config(
    page_title="SSF Calculation",
    page_icon="🚀"
)

# initialization
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done=False
    st.session_state.model=None
# text the title
st.title("SSF Adjustment of FEMA P-695")
# text the introduction
st.header("Introduction")
st.write("Created by Zhiyang Mao and Christianos Burlotos, Stanford, 6/5/2025")
st.write("Supervised by Prof. Gregory Deierlein, Stanford.")
st.write("Last edited: 6/5/2025")
st.write("This website can calculate the spectral shape factor (SSF) using epsilon and Sa Ratio, according to FEMA P-695. User can create the plot based on the input parameters and download files from the analysis.")
st.warning("""
⚠️ This tool is for academic and research purposes only.
Results are for reference only and should not be used for engineering design,
commercial decisions, or risk analysis.
The author assumes no responsibility for any consequences arising from its use.
""")
# FLAG to make sure the input is correct
FLAG_input=1
FLAG_analysis=1

# text the input para
# title
st.header("Input Parameters")
# Fundamental Period
T_input=st.text_input("T",max_chars=100,help="Fundamental Period of the Strucutre")
if T_input != "":
    try:
       T=float(T_input)
       FLAG_input=1
    except ValueError: 
        st.error("Please enter a number")
        FLAG_input=0
else:
   FLAG_input=0
if st.session_state.model!=None and T !=st.session_state.model.T:
    FLAG_analysis=0
# Longitude
Longitude_input=st.text_input("Longitude",max_chars=100,help="[ -125, -65 ]")
if Longitude_input!="":
    try:
       Longitude=float(Longitude_input)
       FLAG_input=1
    except ValueError: 
        st.error("Please enter a number")
        FLAG_input=0
    if Longitude<-125 or Longitude>-65:
        st.error("Longitude out of range. Allowable range: [ -125, -65 ]")
        FLAG_input=0
else:
   FLAG_input=0 
if st.session_state.model!=None and Longitude !=st.session_state.model.Longitude:
   FLAG_analysis=0
# Latitdue
Latitude_input=st.text_input("Latitude", max_chars=100,help=" [ 24.4, 50 ]")
if Latitude_input!="":
    try:
       Latitude=float(Latitude_input)
       FLAG_input=1
    except ValueError: 
        st.error("Please enter a number")
        FLAG_input=0
    if Latitude<24.4 or Latitude>50:
        st.error("Latitude out of range. Allowable range: [ 24.4, 50 ]")
        FLAG_input=0
else:
   FLAG_input=0 
if st.session_state.model!=None and Latitude !=st.session_state.model.Latitude:
    FLAG_analysis=0
# Vs 30
Vs30_input=st.text_input("Vs30 (m/s)",value=260)
if Vs30_input!="":
    try:
       Vs30=float(Vs30_input)
       FLAG_input=1
    except ValueError: 
        st.error("Please enter a number")
        FLAG_input=0
else:
   FLAG_input=0 
if st.session_state.model!=None and Vs30 !=st.session_state.model.Vs30:
    FLAG_analysis=0

# Mechanism
Mechanism_projection={
    "Strike Slip":"SS",
    "Normal Slip":"NS",
    "Reverse slip":"RS",
    "Unspecified":"U"
}
Selected_Mechanism=st.selectbox(
    label="Mechansim of the Main Fault",
    options=("Strike Slip","Normal Slip","Reverse slip","Unspecified"),
    index=3
    )
Mechanism=Mechanism_projection[Selected_Mechanism]
if st.session_state.model!=None and Mechanism !=st.session_state.model.Mechanism:
    FLAG_analysis=0


# Region
Region_Projection={
    "Global":"global",
    "California":"california",
    "China":"china",
    "Italy":"italy",
    "Japan":"japan",
    "New Zealand":"new_zealand",
    "Taiwan":"taiwan",
    "Turkey":"turkey"
}
Selected_Region=st.selectbox(
    label="Region of the Site",
    options=("Global","California","China","Italy","Japan","New Zealand","Taiwan","Turkey"),
    index=1
    )
Region=Region_Projection[Selected_Region]
if st.session_state.model!=None and Region !=st.session_state.model.Region:
    FLAG_analysis=0

# Record Set
Record_Set=st.selectbox(
    label="Record Set",
    options=("Near Field","Far Field","Combined")
)
if st.session_state.model!=None and Record_Set !=st.session_state.model.Record_Set:
    FLAG_analysis=0
# a and b
a_input=st.text_input(label="a",max_chars=100,help="Lower bond for the Sa Ratio: Ta/T")
if a_input!="":
    try:
       a=float(a_input)
       FLAG_input=1
    except ValueError: 
        st.error("Please enter a number")
        FLAG_input=0
else:
   FLAG_input=0 
if st.session_state.model!=None and a !=st.session_state.model.a:
    FLAG_analysis=0
b_input=st.text_input(label="b",max_chars=100,help="Upper bond for the Sa Ratio: Tb/T")
if b_input!="":
    try:
        b=float(b_input)
        FLAG_input=1
    except ValueError: 
        st.error("Please enter a number")
        FLAG_input=0
else:
    FLAG_input=0 
if st.session_state.model!=None and b !=st.session_state.model.b:
    FLAG_analysis=0

# input Sac Form
Sac_File=st.file_uploader(label="Sa collapse",type=["xlsx"],accept_multiple_files=False)
if Sac_File is not None:
    Sac_Form=pd.read_excel(Sac_File)
    st.success("File uploaded successfully")
else:
    st.info("Please upload the file")
    FLAG_input=0
if FLAG_input==0:
    st.stop()
if st.session_state.model!=None and not Sac_Form.equals(st.session_state.model.SacForm):
    FLAG_analysis=0
st.header("Run Analysis")


if st.button("Run"):
    with st.spinner("Analysis Running..."):
        model=SSFAdjustment(T,Longitude,Latitude,Vs30,Mechanism,Region,Record_Set,a,b,Sac_Form)
        Result=model.Run()
        st.session_state.model=model
        st.session_state.analysis_done=True
    st.write("Analysis Complete!")
    FLAG_analysis=1
if FLAG_analysis==0:
    st.stop()

# Plots
if st.session_state.analysis_done == False:
    st.stop()
if "plot_button" not in st.session_state:
    st.session_state.plot_button=None
st.session_state.plot_button=0
    
st.header("Plots")
st.write("What plots do you want to create?")
option_Sac_epsilon=st.checkbox(r"$Sa_{collapse}$-$\epsilon$ Relationship")
option_Sac_Sa_Ratio=st.checkbox(r"$Sa_{collapse}$-$Sa\ Ratio$ Relationship")
option_Fragility_curve=st.checkbox(r"Fragility Curves")
if option_Fragility_curve:
    option_fc_type=st.selectbox(
        label="Which type of fragility curve do you want?",
        options=("epsilon shift","Sa Ratio shift","both"),
        index=2
    )
option_Sac=st.checkbox(r"Spectrum Acceleration at $Sa_{collpase}$")
if st.button("Plot"):
    st.session_state.plot_button=1
if st.session_state.plot_button:
    if option_Sac_epsilon:
        fig=st.session_state.model.SacEpsilonPlot()
        st.pyplot(fig)
    if option_Sac_Sa_Ratio:
        fig=st.session_state.model.SacSaRatioPlot()
        st.pyplot(fig)
    if option_Fragility_curve:
        if option_fc_type=="epsilon shift":
            fig=st.session_state.model.FragilityCurve("Epsilon")
        elif option_fc_type=="Sa Ratio shift":
            fig=st.session_state.model.FragilityCurve("Sa Ratio")
        else:
            fig=st.session_state.model.FragilityCurve("both") 
        st.pyplot(fig)
    if option_Sac:
        fig=st.session_state.model.SacPlot()
        st.pyplot(fig)

st.header("Summary")
st.write(r"SSF for $\epsilon$: {}".format(round(st.session_state.model.SSF_epsilon,3)))
st.write(r"SSF for $Sa\ Ratio$: {}".format(round(st.session_state.model.SSF_Sa_Ratio,3)))
st.header("Files Downlaod")
st.download_button("UHS data",st.session_state.model.UHS.to_csv(index=False),"UHS_data.csv","txt/csv")
st.download_button("epsilon for each record",st.session_state.model.epsilon_for_record.to_csv(index=False),"epsilon_for_each_record.csv","txt/csv")
st.download_button("Sa Ratio for each record",st.session_state.model.Sa_Ratio_for_record.to_csv(index=False),"Sa_Ratio_for_each_record.csv","txt/csv")
st.download_button("Regression of Sa collapse and epsilon",st.session_state.model.epsilon_Regression_Result.to_csv(index=False),"Regression_of_Sa_collapse_and_epsilon.csv","txt/csv")
st.download_button("Regression of Sa collapse and Sa Ratio",st.session_state.model.epsilon_Regression_Result.to_csv(index=False),"Regression_of_Sa_collapse_and_epsilon.csv","txt/csv")




