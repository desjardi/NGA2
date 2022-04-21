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
import re
import os
pd.options.plotting.backend = "plotly"

# Process here the input file to obtain simulation parameters
with open('input') as fobj:
    for line in fobj:
        line_data = re.split(':',line)
        if line_data[0].rstrip()=='Liquid density':
            rho_l=float(line_data[1])
        if line_data[0].rstrip()=='Liquid dynamic viscosity':
            mu_l=float(line_data[1])
        if line_data[0].rstrip()=='Surface tension coefficient':
            sigma=float(line_data[1])
        if line_data[0].rstrip()=='Droplet radius':
            Rinit=float(line_data[1])
        if line_data[0].rstrip()=='Static contact angle':
            CA0=float(line_data[1])
        if line_data[0].rstrip()=='Hole size':
            hs=float(line_data[1])
        if line_data[0].rstrip()=='Hole dist':
            hd=float(line_data[1])
# Compute Rayleigh time scale
Tref=(rho_l*Rinit**3/sigma)**(1/2)
# Compute porosity
eps=hs**2/hd**2
# Dry Cassie-Baxter angle
DCB=(1-eps)*math.cos(CA0)-eps
# Wet Cassie-Baxter angle
WCB=(1-eps)*math.cos(CA0)+eps

# Define imbibed volume figure
df = pd.read_csv('monitor/dropinfo', delim_whitespace=True, header=None, skiprows=2, usecols=[1, 4, 5], names=['Time', 'Vtot', 'Vimb'])
df['Iper'] = 100*df['Vimb']/df['Vtot']
df['Tper'] = 100*df['Vtot']/df['Vtot']
df['NTime'] = df['Time']/Tref

fig1=go.Figure()
fig1.add_trace(go.Scatter(x=df['NTime'], y=df['Iper'], fill='tozeroy', mode='none', showlegend=False)) # fill down to xaxis
fig1.add_trace(go.Scatter(x=df['NTime'], y=df['Tper'], fill='tonexty', mode='none', showlegend=False)) # fill to trace0 y
fig1.update_layout(width=800,height=600)
#fig1.update_layout(title_text='Imbibition over time',title_font_size=36,title_x=0.5)
fig1.update_xaxes(title_text='Normalized time',title_font_size=24,tickfont_size=24)
fig1.update_yaxes(title_text='Percent imbibed',title_font_size=24,tickfont_size=24,range=[0,100])
fig1.add_annotation(x=1, y=90,text='Non-imbibed',showarrow=False, font_size=24, font_color='darkslategrey')
fig1.add_annotation(x=4, y=10,text='Imbibed'    ,showarrow=False, font_size=24, font_color='darkslategrey')

# Prepare radius=f(height) here
rf = pd.read_csv('radius/radius_0.00000E+00', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['y'])
rf0 = rf.copy()
rtime=[]
for radfile in os.listdir('radius'):
    rad_name = re.split('_',radfile)
    rtime.append(float(rad_name[1]))
rtime.sort()
ftime=[]
for t in rtime:
    formatted_t = "{:12.5E}".format(t)
    ftime.append(formatted_t.lstrip())
for ft in ftime:
    radfile='radius/radius_'+ft
    tmp = pd.read_csv(radfile, delim_whitespace=True, header=None, skiprows=1, usecols=[1], names=[ft])
    rf[ft] = tmp[ft]

# Now generate the wetted radius and contact angle
wet_rad=[]
angle=[]
cosangle=[]
Cap=[]
for myind in range(len(rtime)):
    myrf=rf0.copy()
    myrf['y']=rf['y']
    myrf['R']=rf[ftime[myind]]
    myrf=myrf[myrf['y']>=0]
    myrf=myrf[myrf['R']>0]
    #print(myind,len(myrf.index))
    # Get wetting radius
    #wet_rad.append(myrf.iat[0,1])  # This is a first order estimate
    try:
        wet_rad.append(1.5*myrf.iat[0,1]-0.5*myrf.iat[1,1])  # This is a second-order estimate
        # Get contact angle
        myangle=math.atan2(myrf.iat[1,1]-myrf.iat[0,1],myrf.iat[1,0]-myrf.iat[0,0])+0.5*math.pi
        angle.append(myangle)
        cosangle.append(math.cos(myangle))
        # Get CL velocity
        if myind==0:
            vel=0
        else:
            vel=(wet_rad[myind]-wet_rad[myind-1])/(rtime[myind]-rtime[myind-1])
        myCap=mu_l*vel/sigma
        Cap.append(myCap)
    except Exception as e:
        break

wr = pd.DataFrame(list(zip(rtime,wet_rad,angle,cosangle,Cap)),columns=['Time','R','angle','cosangle','Cap'])
fig3=go.Figure()
fig3.add_trace(go.Scatter(x=wr['Time']/Tref, y=wr['R']/Rinit, mode='lines', showlegend=False))
fig3.update_layout(width=800,height=600)
fig3.update_xaxes(title_text='Normalized time',type='log',title_font_size=24,tickfont_size=24,range=[math.log10(0.01),math.log10(2)])
fig3.update_xaxes(tickvals=[0.0001,0.001,0.01,0.1,1,10,100])
fig3.update_yaxes(title_text='Normalized wetted radius',type='log',title_font_size=24,tickfont_size=24,range=[math.log10(0.07),math.log10(2)])
fig3.update_yaxes(tickvals=[0.0001,0.001,0.01,0.1,1,10,100])

fig4=go.Figure()
fig4.add_trace(go.Scatter(x=wr['Time']/Tref, y=wr['cosangle'], mode='lines', showlegend=False))
fig4.update_layout(width=800,height=600)
fig4.update_xaxes(title_text='Normalized time',title_font_size=24,tickfont_size=24)
fig4.update_yaxes(title_text='Cosine of contact angle',title_font_size=24,tickfont_size=24,range=[-1,+1])
fig4.add_shape(type='line',x0=0,y0=math.cos(CA0/180*math.pi),x1=rtime[len(rtime)-1]/Tref,y1=math.cos(CA0/180*math.pi),line_color='red')
fig4.add_annotation(x=0.7,y=math.cos(CA0/180*math.pi)+0.1,text='Static contact',showarrow=False,font_size=16,font_color='red')

fig5=go.Figure()
fig5.add_trace(go.Scatter(x=wr['Cap'], y=wr['cosangle'],mode='lines+markers',showlegend=False,line=dict(color='black',dash='dot',width=1),marker=dict(color=wr['Time']/Tref,colorscale='Viridis',line_width=1)))
fig5.update_layout(width=800,height=600)
fig5.update_xaxes(title_text='Capillary number',title_font_size=24,tickfont_size=24)
fig5.update_yaxes(title_text='Cosine of contact angle',title_font_size=24,tickfont_size=24,range=[-1,+1])
fig5.add_shape(type='line',x0=min(Cap),y0=math.cos(CA0/180*math.pi),x1=max(Cap),y1=math.cos(CA0/180*math.pi),line_color='darkslategrey')
#fig5.add_annotation(x=0.5*min(Cap),y=math.cos(CA0/180*math.pi)+0.1,text='Static contact',showarrow=False,font_size=16,font_color='darkslategrey')
fig5.add_shape(type='line',x0=min(Cap),y0=math.cos(WCB),x1=max(Cap),y1=math.cos(WCB),line_color='blue')
fig5.add_shape(type='line',x0=min(Cap),y0=math.cos(DCB),x1=max(Cap),y1=math.cos(DCB),line_color='red')

# This is where we define the dashboard layout
app = dash.Dash(__name__)
app.layout = html.Div(style={"margin-left": "15px"},children=[

    # Title of doc
    dcc.Markdown('''# Imbibition Project'''),
    dcc.Markdown('''*NGA2 Dashboard written by O. Desjardins, last updated 01/18/2021*'''),
    
    # Intro
    dcc.Markdown('''
    ## Overview

    In this dashboard, we post-process the raw data generated by NGA2's imbibition
    case. This simulation is based on the Sahoo and Louge experiment of droplet imbibition
    and spreading on a perforated plate, conducted in late 2018 on the ISS.
    '''),
    
    # Parameters
    dcc.Markdown(f'''
    ---
    ## Simulation parameters

    By analyzing the input file, we have detected the following parameters:
    - Surface tension coefficient: \u03C3 = {sigma}
    - Droplet radius: R = {Rinit}
    - Static contact angle: CA = {CA0}
    - Porosity of the plate: \u03B5 = {eps}
    
    The resulting Rayleigh time scale is T = {Tref}
    '''),
    
    # Imbibed volume over time
    dcc.Markdown(f'''
    ---
    ## Imbibed liquid volume over time
    
    The graph below shows the fraction of the droplet volume imbibed in the
    perforated plate over time. Time has been normalized by the Rayleigh timescale.
    '''),
    dcc.Graph(id='Volume_imbibed',figure=fig1),
    
    # Radius=f(y,t)
    dcc.Markdown(f'''
    ---
    ## Interface shape over time
    
    The graph below shows the orthoradially-averaged interface shape as a function of time.
    Time has been normalized by the Rayleigh timescale.
    '''),
    dcc.Graph(id='Radius-slider'),
    daq.Slider(id='time-slider',min=0,max=len(rtime)-1,step=1,value=0,updatemode='drag',size=700,marks={
        0: {'label': str(rtime[0]), 'style': {'color': 'black'}},
        len(rtime)-1: {'label': str(rtime[len(rtime)-1]), 'style': {'color': 'black'}}
    }),
    html.Div(id='time-slider-output-container'),
    html.Br(),
    
    # Wet radius over time
    dcc.Markdown(f'''
    ---
    ## Wetted liquid radius over time
    
    The graph below shows the normalized wetted liquid radius as a function of normalized time.
    '''),
    dcc.Graph(id='Wetted_radius',figure=fig3),
    
    # Contact angle over time
    dcc.Markdown(f'''
    ---
    ## Contact angle over time
    
    The graph below shows the cosine of the contact angle as a function of normalized time.
    '''),
    dcc.Graph(id='Contact_angle',figure=fig4),
    
    # Contact angle / Ca number trajectory
    dcc.Markdown(f'''
    ---
    ## Trajectory in capillary number/contact angle space
    
    The graph below shows the trajectory of the simulation in capillary number/contact angle space.
    '''),
    dcc.Graph(id='Contact_Capillary',figure=fig5),
    
])

@app.callback(
    dash.dependencies.Output('Radius-slider','figure'),
    [dash.dependencies.Input('time-slider','value')])
def update_figure(value):
    myind=value
    myrf=rf0.copy()
    myrf['y']=rf['y']/Rinit
    myrf['R']=rf[ftime[myind]]/Rinit
    myrf=myrf[myrf['y']>=0]
    myrf=myrf[myrf['R']>0]

    fig2=go.Figure()
    fig2.add_trace(go.Scatter(x=myrf['R'], y=myrf['y'], showlegend=False, mode='lines+markers', line_color='red', line_shape='spline'))
    fig2.update_layout(width=700,height=800)
    fig2.update_xaxes(title_text='Normalized distance from centerline axis of drop',title_font_size=24,tickfont_size=24,range=[0,2.5])
    fig2.update_yaxes(title_text='Normalized height above plate',title_font_size=24,tickfont_size=24,range=[0,3])
    fig2.add_annotation(x=1.5, y=2.75,text='Normalized time={:.3f}'.format(rtime[myind]/Tref),showarrow=False, font_size=24, font_color='darkslategrey')

    return fig2


if __name__ == '__main__':
    app.run_server(debug=True)
