# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_daq as daq
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import math
import numpy as np
import re
import os
pd.options.plotting.backend = "plotly"

# Create the dashboard
app = dash.Dash(__name__)


def get_cond(rho_,T_):
    A1=-105.161
    A2=+0.9007
    A3=+0.0007
    A4=+3.50e-15
    A5=+3.76e-10
    A6=+0.7500
    A7=+0.0017
    # Evaluate conductivity in mW/(m.K)
    cond=(A1+A2*rho_+A3*np.power(rho_,2)+A4*np.power(rho_,3)*np.power(T_,3)+A5*np.power(rho_,4)+A6*T_+A7*np.power(T_,2))/np.power(T_,0.5)
    k_=0.001*cond
    return k_


# Define here some parameters
Rcst=8.314           # J/(mol.K)
Wmlr=44.01e-3        # kg/mol
Cp=40.0/Wmlr         # J/(kg.K)
Cv=Cp-Rcst/Wmlr      # J/(kg.K)
Gamma=Cp/Cv          # -
Cs=500.0             # J/(kg.K)
k_steel=45           # W/(kg.m)
Tinlet=430           # K
Tinit =300           # K
MFR=0.2              # kg/m^3
Minit=9.87217        # kg
Tw=300               # K


# Vessel geometry
radius=0.4           # m
length=2.5           # m
area=2.0*math.pi*radius*length+2.0*math.pi*radius*radius
Vtotal=math.pi*radius*radius*length # m^3
m_steel=268.5        # kg of steel inside
L_steel=0.0762       # m of steel thickness
dx=2.6/80

# Temperature graph
def create_Tfig():
    
    # Create fig
    Tfig=go.Figure()
    
    # Read and plot the Farther Farms data
    df=pd.read_csv('FFdata/temperature.txt',delim_whitespace=True,header=None,skiprows=0,usecols=[0,1],names=['Time','Temp'])
    Tfig.add_trace(go.Scatter(name='FF experiment - with sled',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(width=4)))
    
    # Some formatting
    Tfig.update_layout(width=1200,height=800)
    Tfig.update_xaxes(title_text='Time (min)',title_font_size=24,tickfont_size=24,range=[0,20])
    Tfig.update_yaxes(title_text='Temperature (K)',title_font_size=24,tickfont_size=24)#,range=[280,450])
    #Tfig.add_shape(type='line',x0=0,y0=Tinit,x1=df['Time'].iloc[-1]/60,y1=Tinit,line_color='black')
    #Tfig.add_annotation(x=7.5,y=Tinit-7,text='Tinit',showarrow=False,font_size=16,font_color='black')
    #Tfig.add_shape(type='line',x0=0,y0=Tinlet,x1=df['Time'].iloc[-1]/60,y1=Tinlet,line_color='green')
    #Tfig.add_annotation(x=7.5,y=Tinlet+7,text='Tinlet',showarrow=False,font_size=16,font_color='green')
    Tfig.update_layout(legend=dict(font=dict(size=14)))
    
    # Add temperature from monitor_basket_and_wallmodel
    df=pd.read_csv('monitor_basket_and_wallmodel/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5,6],names=['Time','Temp','Mass','Pres','Twall'])
    Tfig.add_trace(go.Scatter(name='No product - Tvessel',x=df['Time']/60,y=df['Temp'] ,mode='lines',showlegend=True,line=dict(color='darkgreen',width=2)))
    Tfig.add_trace(go.Scatter(name='No product - Twall',  x=df['Time']/60,y=df['Twall'],mode='lines',showlegend=True,line=dict(color='red',width=2,dash='dot')))
    
    # Add temperature running now
    df=pd.read_csv('monitor/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5,6],names=['Time','Temp','Mass','Pres','Twall'])
    Tfig.add_trace(go.Scatter(name='With product - Tvessel',x=df['Time']/60,y=df['Temp'] ,mode='lines',showlegend=True,line=dict(color='black',width=2)))
    Tfig.add_trace(go.Scatter(name='With product - Twall',  x=df['Time']/60,y=df['Twall'],mode='lines',showlegend=True,line=dict(color='black',width=2,dash='dot')))
    
    
    # Read and plot the Farther Farms data
    #df=pd.read_csv('FFdata/temp_nosled.txt',delim_whitespace=True,header=None,skiprows=0,usecols=[0,1],names=['Time','Temp'])
    #Tfig.add_trace(go.Scatter(name='FF experiment - no sled',x=(df['Time']-40)/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(width=4)))
    
    
    # 0D model of the vessel using the measured Tinlet
    df=pd.read_csv('FFdata/inlet_temp.txt',delim_whitespace=True,header=None,skiprows=0,usecols=[0,1],names=['Time','Tin'])
    
    myTime=df['Time'].tolist()
    myTin=df['Tin'].tolist()
    myMass=np.zeros(len(myTime))
    myMass[0]=Minit
    myTemp=np.zeros(len(myTime))
    myTemp[0]=Tinit
    myTemp2=np.zeros(len(myTime))
    myTemp2[0]=Tinit
    
    for n in range(0,len(myTime)-1):
        myMass[n+1]=myMass[n]+(myTime[n+1]-myTime[n])*MFR
        myTemp[n+1]=((myMass[n]+m_steel*Cs/Cv)*myTemp[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*Tinlet))/(myMass[n+1]+m_steel*Cs/Cv)
        #\myTemp2[n+1]=((myMass[n]+m_steel*Cs/Cv)*myTemp2[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*myTin[n]-2*area/(Cv*L_steel)*(myTemp2[n]-Tw)))/(myMass[n+1]+m_steel*Cs/Cv)
        myTemp2[n+1]=((myMass[n]+m_steel*Cs/Cv)*myTemp2[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*myTin[n]-area/(Cv*L_steel)*(myTemp2[n]-Tw)))/(myMass[n+1]+m_steel*Cs/Cv)
        
    
    df=pd.DataFrame(list(zip(myTime,myTin,myTemp,myTemp2)),columns=['Time','Tin','Temp','Temp2'])
    Tfig.add_trace(go.Scatter(name='Adiabatic model',x=df['Time']/60,y=df['Temp'] ,mode='lines',showlegend=True,line=dict(color='red',width=2)))
    #Tfig.add_trace(go.Scatter(name='Non-adiabatic model',x=df['Time']/60,y=df['Temp2'],mode='lines',showlegend=True,line=dict(color='blue',width=2)))
    
    return Tfig




# Pressure graph
def create_Pfig():
    
    # Create fig
    Pfig=go.Figure()
    
    # Read and plot the Farther Farms data
    df=pd.read_csv('FFdata/pressure.txt',delim_whitespace=True,header=None,skiprows=0,usecols=[0,2],names=['Time','Pres'])
    Pfig.add_trace(go.Scatter(name='Farther Farms experiment',x=df['Time'],y=df['Pres'],mode='markers',showlegend=True,line=dict(width=2)))
    
    Pfig.update_layout(width=1200,height=800)
    Pfig.update_xaxes(title_text='Time (min)',title_font_size=24,tickfont_size=24,range=[0,20])
    Pfig.update_yaxes(title_text='Pressure (bar)',title_font_size=24,tickfont_size=24)
    #Pfig.add_shape(type='line',x0=0,y0=Tinit,x1=df['Time'].iloc[-1]/60,y1=Tinit,line_color='black')
    #Pfig.add_annotation(x=7.5,y=Tinit-7,text='Tinit',showarrow=False,font_size=16,font_color='black')
    #Pfig.add_shape(type='line',x0=0,y0=Tinlet,x1=df['Time'].iloc[-1]/60,y1=Tinlet,line_color='green')
    #Pfig.add_annotation(x=7.5,y=Tinlet+13,text='Tinlet',showarrow=False,font_size=16,font_color='green')
    Pfig.update_layout(legend=dict(font=dict(size=14)))
    
    # Read and plot current data
    df=pd.read_csv('monitor/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5,6],names=['Time','Temp','Mass','Pres','Twall'])
    Pfig.add_trace(go.Scatter(name='NGA2 simulation - 1 min shift',x=df['Time']/60+1,y=df['Pres']*1e-5,mode='lines',showlegend=True,line=dict(width=2)))
    

    return Pfig




# This is where we define the dashboard layout
def serve_layout():
    return html.Div(style={"margin-left": "15px"},children=[
    
    # Title of doc
    dcc.Markdown('''# Farther Farms Project'''),
    dcc.Markdown('''*NGA2 Dashboard written by O. Desjardins, last updated 02/06/2021*'''),
    
    # Intro
    dcc.Markdown('''
    ## Overview
    In this dashboard, we post-process the raw data generated by NGA2's pvessel
    case. This simulation is based on an experiment done by Farther Farms where
    a pressure vessel is filled with heated CO2.
    '''),
    
    # Temperature over time
    dcc.Markdown(f'''
    ---
    ## Average temperature in the vessel
    The graph below shows the evolution of the average temperature inside the pressurized vessel.
    '''),
    dcc.Graph(id='Tgraph',figure=create_Tfig()),
    
    # Pressure over time
    dcc.Markdown(f'''
    ---
    ## Pressure in the vessel
    The graph below shows the evolution of the pressure inside the vessel.
    '''),
    dcc.Graph(id='Pgraph',figure=create_Pfig()),
])


# This is where we set the layout and run the server
app.layout = serve_layout
if __name__ == '__main__':
    app.run_server(debug=True)
