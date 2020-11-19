# # ToDos

# ToTal Costs überprüfen
# Auswertung mit und ohne WP
# Auswertung mit und ohne untere Caverne
# Mit kleinerer Caverne und ohne MinCapacities
# Reducted Costs?
# Nur Pumped Hydro
# Anderer Fernwärmepreis?

#Erledigt:
# Abwärmewert ausgeben und validieren
# Wärmepumpe: Wärmebilanz: - Verluste + Strom - Entnahmen in die Fernwärme


using JuMP, XLSX, CPLEX

#m = Model(CPLEX.Optimizer) hier kann man keine CPLEX routinen aufrufen!
m = direct_model(CPLEX.Optimizer()) # Nur so geht es scheinbar
mm = backend(m)


set_optimizer_attribute(m, "CPXPARAM_WorkDir", "C:\\Gerhard\\Model\\Julia\\SeasonalStorage")
CPXsetlogfilename( mm.env, "mylogfile.txt", "a")# a for append 


set_optimizer_attribute(m, "CPXPARAM_ScreenOutput", 1)   # Enable output
set_optimizer_attribute(m, "CPXPARAM_SolutionType", 2)  # No crossover
set_optimizer_attribute(m, "CPXPARAM_LPMethod", 4)      # Use barrier

#set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)  # No presolve
#set_optimizer_attribute(m, "CPXPARAM_TimeLimit", 3600)      # Time limit in seconds
#set_optimizer_attribute(m, "CPXPARAM_Threads", 8)      # Upper limit for the number of threads
#set_optimizer_attribute(m, "CPXPARAM_LPMethod", 4)      # Use barrier

"""
    cop(tC::Union{Int64, Float64}, tH::Union{Int64, Float64}, eff::Float64)

DOCSTRING

# Arguments:
- `tC`: DESCRIPTION
- `tH`: DESCRIPTION
- `eff`: DESCRIPTION
"""
function cop(tC::Union{Int64,Float64},tH::Union{Int64,Float64},eff::Float64)
    # We assume that the heat source tC is only cooled a little = 5°C and tH is the target temperature
    # eff: is the faktor of real life efficiency to the heat pumpt carnot efficiency 
    # the carnot efficiency is calculated in Kelvin!!
    # https://de.wikipedia.org/wiki/Carnot-Wirkungsgrad#Analoge_Gr%C3%B6%C3%9Fen_f%C3%BCr_W%C3%A4rmepumpen_und_K%C3%A4ltemaschinen
    eff*(tH+273.15)/(tH-tC)
end

"""
    annualizedCost_EuropMWpYr(invCost_EuroPkW_elec::Float64, depreciationYears::Int64, wacc::Float64)

DOCSTRING

# Arguments:
- `invCost_EuroPkW_elec`: DESCRIPTION
- `depreciationYears`: DESCRIPTION
- `wacc`: DESCRIPTION
"""
function annualizedCost_EuropMWpYr(invCost_EuroPkW_elec::Float64,depreciationYears::Int64, wacc::Float64)
    #wacc in percent: [0-100]
    #https://de.wikipedia.org/wiki/Rentenrechnung
    # Barwert Nachschüssig: JahresRate=Bnach*q^n*i/(q^n-1)
    wacc/100*(1+wacc/100)^depreciationYears/((1+wacc/100)^depreciationYears-1)*invCost_EuroPkW_elec*1000
end

"""
    t_m (delta_t_leftside, delta_t_rightside)
    effective mean temperature difference between cold and hot side of a heat exchanger 
    delta_t_leftside=temperature differenc between the two water streams at the left side 
    delta_t_rightside=temperature differenc between the two water streams at the right side 
    formula works for both gegenstom und  gleichstrom wärmetauscher
    https://de.wikipedia.org/wiki/W%C3%A4rmetauscher#Berechnung_und_Bewertung_von_Rekuperatoren
"""
function t_m(delta_t_leftside::Union{Int64, Float64}, delta_t_rightside::Union{Int64, Float64})
    if delta_t_leftside!=delta_t_rightside 
        (delta_t_leftside-delta_t_rightside)/log(delta_t_leftside/delta_t_rightside) 
    else delta_t_leftside end
end

log

"Binary_period: is a time range where the use of caverns is correctly modelled with binary variables"
function add_cavern_var!(m::Model,cav_levels::Tuple{String,String},temp_levels::OrdinalRange{Int64,Int64},VL::Array{Float64,1},
    CapLimitFW_Waermetauscher_MW::Float64,
    CapLimitWasteH_Waermetauscher_MW::Float64,
    CapLimitRiver_Waermetauscher_MW::Float64,
    riverTemp::Array{Float64,1},
    time_steps::OrdinalRange{Int64,Int64})



    # Durch die NumberOfCavernsPerLevel wird die gesamte Speichergröße optimiert. Wir nehmen eine fixe Größe für eine Röhre an (capOneCavern_1000m3).
    # @variable(m, NumberOfCavernsPerLevelC[lv in cav_levels],lower_bound=(if lv=="hi" 30 else 0 end), upper_bound=200);
    # @variable(m, 200>=NumberOfCavernsPerLevelD[lv in cav_levels]>=(if lv=="hi" 30 else 0 end));
    @variable(m, (if lv=="hi"; 10 else 10 end)<=NumberOfCavernsPerLevel[lv in cav_levels]<=200);#(if lv=="hi"; 60 else  end)


    # Number of Caverns needed for each temperature level and cavern level at the given timestep. It is assumed that in each single cavern 
    @variable(m, 0<=NrOfCavernNeedAtTempLNow[cav_levels,temp_levels,binary_period;length(binary_period)>0]<=200, integer=true);

    # M3_storage is the storage content at the begin of the time_steps
    @variable(m, M3_storage[cav_levels,temp_levels,time_steps]>=0);

    # m3 direkt Einspeisung aus dem Saisonalspeicher in die Fernwärme geht nur wenn VL[t]+temp_levels.step=tl dem Speicherlevel entspricht 
    # die VL Tem wird auf die temp_levls gerundet
    # wenn VL[t]+temp_levels.step>temp_levels.stop dann wird bei 95° entnommen und mit der Wärmepumpe nachgeheitz bis zur VL[t]+temp_levels.step

    @variable(m, 0.2<=heatExchangerCap_CavToFW_MW<=CapLimitFW_Waermetauscher_MW);
    @variable(m, 0.2<=heatExchangerCap_Wasteheat_MW<=CapLimitWasteH_Waermetauscher_MW);
    @variable(m, 0.2<=heatExchangerCap_river_MW<=CapLimitRiver_Waermetauscher_MW);

    @variable(m, 0<=districtHeatFromStorage_m3pS[tl in temp_levels,t in time_steps;VL[t]+temp_levels.step<=tl]);

    @variable(m, 0<=riverwater_m3pS[t in time_steps;riverTemp[t]>=(temp_levels.start+2*temp_levels.step)]);
    @variable(m, 0<=returnRiverwater_m3pS[tl in temp_levels,t in time_steps;tl+2*temp_levels.step<=riverTemp[t]]);

    #=
    NumberOfCavernsPerLevel=m[:NumberOfCavernsPerLevel]
    NrOfCavernNeedAtTempLNow=m[:NrOfCavernNeedAtTempLNow]
    M3_storage=m[:M3_storage]
    districtHeatFromStorage_m3pS=m[:districtHeatFromStorage_m3pS]
    =m[:]
    =#
    return m
end


#binary_period is a time period where caverns are simulated with integers
function add_cavern_conUexp!(m::Model,cav_levels::Tuple{String,String},capOneCavern_1000m3::Int64,
    invCostCavern_EuroPm3::Float64,depreciationYearsCavern::Int64, 
    costHeatExchanger_wasteheat_EuroPkWth::Float64,depreciationYearsHeatEx::Int64,
    wacc::Float64, 
    temp_levels::OrdinalRange{Int64,Int64},time_steps::OrdinalRange{Int64,Int64},binary_period::OrdinalRange{Int64,Int64})

    # Variablen Definitionen die gebraucht werden aus anderen Funktionen laden
    NumberOfCavernsPerLevel=m[:NumberOfCavernsPerLevel]
    if length(binary_period)>0; NrOfCavernNeedAtTempLNow=m[:NrOfCavernNeedAtTempLNow] end
    M3_storage=m[:M3_storage]
    heatExchangerCap_CavToFW_MW=m[:heatExchangerCap_CavToFW_MW]
    heatExchangerCap_Wasteheat_MW=m[:heatExchangerCap_Wasteheat_MW]
    heatExchangerCap_river_MW=m[:heatExchangerCap_river_MW]
    returnRiverwater_m3pS=m[:returnRiverwater_m3pS]
    riverwater_m3pS=m[:riverwater_m3pS]
 
    @constraint(m, storageCap[h in cav_levels,t in time_steps], sum(M3_storage[h,:,t])<=NumberOfCavernsPerLevel[h]*capOneCavern_1000m3*1000);
    @constraint(m, cavernRequiredLimit[h in cav_levels,t in binary_period; length(binary_period)>0 ], sum(NrOfCavernNeedAtTempLNow[h,tl,t] for tl in temp_levels)<=NumberOfCavernsPerLevel[h]);
    @constraint(m, cavernUseLimit[h in cav_levels,tl in temp_levels,t in binary_period; length(binary_period)>0 ], M3_storage[h,tl,t]<=NrOfCavernNeedAtTempLNow[h,tl,t]*capOneCavern_1000m3*1000)

    @constraint(m,riverBal[t in time_steps; riverTemp[t]>=(temp_levels.start+2*temp_levels.step)],
        #Warum soll die gleichung stimmen? Durch den Wärmetauscher sind die Massenströme ja entkoppelt, 
        # wenn das delta t im Wärmetauscher aber auf beiden seiten gleich groß ist müssen die Massenströme auch gleich sein
        # Im Speicher selber muss ja sowieso gelten dass das Wasser erhlaten ist. riverwater_m3pS ist eigentlich überflüssig
        sum(returnRiverwater_m3pS[tl,t] for tl in temp_levels if (tl+2*temp_levels.step)<=riverTemp[t])==riverwater_m3pS[t])
    
    @expression(m,totalCostCavern_exp, 
    (
    sum(NumberOfCavernsPerLevel[:])*capOneCavern_1000m3*annualizedCost_EuropMWpYr(invCostCavern_EuroPm3,depreciationYearsCavern,wacc)
    +(heatExchangerCap_CavToFW_MW +heatExchangerCap_Wasteheat_MW+heatExchangerCap_river_MW)*annualizedCost_EuropMWpYr(costHeatExchanger_wasteheat_EuroPkWth,depreciationYearsHeatEx,wacc)
    ))

    #= 
    totalCostCavern_exp=m[:totalCostCavern_exp]
    =#
    return m
end




#spezifische Wärme c von Wasser 4.180 kJ/(kg K)
#Q_MW=m3pS*1000*deltaT*4.180/1000 in MW 
#--> Q_MW=m3pS*deltaT*4.180
# m3pS=Q_MW/deltaT/4.180

# Q_WP_ab=COP*Pelec
# verl=Q_WP_ab/(Qzu+Pelec)   z.B. verl=0.97

# Nachprüfen1:
# verl=Q_WP_ab/(Qzu+Pelec) --> Qzu =Q_WP_ab/verl-Pelec
# Qzu=COP*Pelec/verl-Pelec=(COP/verl-1)*Pelec

# -->
# Pelec= Qzu/(COP/verl-1)	
# Q_WP_ab= Qzu/(1/verl-1/COP)
# Qzu=(COP/verl-1)*Pelec    
# m_zu=Qzu/(c*dT) 
# m_ab=Qab/(c*dT_FW)

## Nachprüfen2:
## Q_WP_ab=Qzu/(1/verl-1/COP)
## Q_WP_ab= (COP/verl-1)*Pelec/(1/verl-1/COP)
## ((COP-verl)/verl)/((COP-verl)/(COP*verl))=(1/verl)/(1/(COP*verl))=COP
## -->Q_WP_ab=COP*Pelec

# Q_WP_zu_MW=m3pS_zu*deltaT_zu*4.180
# Q_WP_ab_MW=m3pS_zu*deltaT_zu*4.180/(1/verl-1/COP)
# m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP)=im Speicher wenn beide delta Ts gleich sind =m3pS_zu*1/(1/verl-1/COP)


# m3pS_ab=HP_m3Ps[tC,tH,t]/(1/verl-1/COP)
# DisctrictHeat_MW=HP_m3Ps[tC,tH,t]*deltaT_zu*4.180/(1/verl-1/COP)
# Pelec_MW=HP_m3Ps[tC,tH,t]*deltaT_zu*4.180/(COP/verl-1)


function add_heatpumpInCavern_var!(m::Model,temp_levels::OrdinalRange{Int64,Int64},supply_levels::OrdinalRange{Int64,Int64},HP_tH_upperbound::Array{Float64,1},
    time_steps::OrdinalRange{Int64,Int64})
    # A supply level of 100 is the flag to provide the required supply temperature of the district heating system
    # HP_m3Ps[tC,tH,t] ist the m3/s mass flow on the source side and takes source heat from tC and cools it to tC-5 and supplys heat at tH by reheating tH-5
    # on the supply side the mass flow is HP_m3Ps[tC,tH,t]/(1/verl-1/COP)
    # if tH=100: take return flow water from district heat and reheat it to supply temperature 
    
    #Variable for installed HP electrical capacity. 
    @variable(m, 1<=HP_Cap_MW) 

    # HP_m3Ps[tC,tH,t] ist the m3/s mass flow on the source side and takes source heat from tC and cools it to tC-5 and supplys heat at tH by reheating tH-5
    # on the supply side the mass flow is HP_m3Ps[tC,tH,t]/(1/verl-1/COP)
    # HP_tH_upperbound::Float64=[(VL[t]+temp_levels.step>temp_levels.stop) ? 100 : temp_levels.stop for t in time_steps ]
    @variable(m, HP_m3Ps[tC in temp_levels,tH in supply_levels,t in time_steps; temp_levels.start<tC<tH<=HP_tH_upperbound[t] ]>=0)

    #=
    HP_Cap_MW=m[:HP_Cap_MW]
    HP_m3Ps=m[:HP_m3Ps]
    =m[:]
    =#

end



function add_heatpumpInCavern_expCon!(
    m::Model,
    temp_levels::OrdinalRange{Int64,Int64},
    supply_levels::OrdinalRange{Int64,Int64},
    time_steps::OrdinalRange{Int64,Int64},
    VL::Array{Float64,1}, 
    HP_tH_upperbound::Array{Float64,1},
    RL::Array{Float64,1}, 
    invCostHP_EuroPkW_elec::Float64,
    depreciationYearsHP::Int64, 
    wacc::Float64
    )


    HP_Cap_MW=m[:HP_Cap_MW]
    HP_m3Ps=m[:HP_m3Ps]


    #für die Ausgabe
    @expression(m,HP_supply_m3Ps_exp[tC in temp_levels, tH in supply_levels, t in time_steps;temp_levels.start<tC<tH<=HP_tH_upperbound[t]],
            (
            #Output m3/s of HP on the supply side
            # m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP) =im Speicher wenn beide delta Ts gleich sind =m3pS_zu*1/(1/verl-1/COP)
            if tH < supply_levels.stop
                HP_m3Ps[tC,tH,t]/(1/verl-1/cop(tC,tH,0.5)) 
            else 0 end 
            + 
            # Nur wenn HP_tH_upperbound[t]==100
            # Wenn tH=100 ist dann wird temp_levels.stop (95°C) auf die VL + 5°C aufgeheizt und dann in einem Wärmetauscher auf RL+5°C abgekühlt und geht wieder in den SPeicher  
            # m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP)
            if tH==supply_levels.stop==HP_tH_upperbound[t]
                HP_m3Ps[tC,tH,t]*temp_levels.step/(VL[t]+temp_levels.step-temp_levels.stop)/(1/verl-1/cop(tC,VL[t]+temp_levels.step,0.5)) 
            else 0 end
            ))
    

    # Storage level balance: here is the heat pump part of the balance, 
    # there is also the pumped hydro part and the gasheater and electric heater part    
    @expression(m,storageLevBal_HP_exp[tl in temp_levels, t in time_steps],
            (
            #Input of the Heatpump on the heat source side: 10<=tl<=95
            #95°C can be cooled to supply at the district heat supply temperature >100 
            if temp_levels.start < tl <= HP_tH_upperbound[t]-temp_levels.step
                -sum(HP_m3Ps[tl,tH,t]*3600 for tH in supply_levels if tl<tH<=HP_tH_upperbound[t])
            else 0 end
            +
            #Output of cooled liquid at the Heatpump heat source side:
            if tl+temp_levels.step<=HP_tH_upperbound[t]-temp_levels.step 
                    sum(HP_m3Ps[tl+temp_levels.step,tH,t]*3600 for tH in supply_levels if tl+temp_levels.step<tH<=HP_tH_upperbound[t])
            else 0 end

            + 
            #Output of HP on the supply side
            # m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP) =im Speicher wenn beide delta Ts gleich sind =m3pS_zu*1/(1/verl-1/COP)
            if supply_levels.start <= tl < supply_levels.stop
                sum(HP_m3Ps[tC,tl,t]*3600/(1/verl-1/cop(tC,tl,0.5)) for tC in temp_levels if temp_levels.start<tC<tl )
            else 0 end 
            + 
            # Nur wenn HP_tH_upperbound[t]=100
            # Wenn tH=100 ist dann wird temp_levels.stop (95°C) auf die VL + 5°C aufgeheizt 
            # und dann in einem Wärmetauscher auf RL+5°C abgekühlt und geht wieder in den SPeicher  
            # m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP)
            if HP_tH_upperbound[t]==supply_levels.stop && tl ==(RL[t]+temp_levels.step)
                sum(HP_m3Ps[tC,supply_levels.stop,t]*3600*temp_levels.step/(VL[t]+temp_levels.step-temp_levels.stop)/(1/verl-1/cop(tC,VL[t]+temp_levels.step,0.5)) for tC in temp_levels if temp_levels.start<tC)
            else 0 end
             
            +
            #Input  of HP on the supply side    
            # m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP) =im Speicher wenn beide delta Ts gleich sind =m3pS_zu*1/(1/verl-1/COP)
            if supply_levels.start <=(tl+temp_levels.step) < supply_levels.stop # >=15°C <100°C : With tH=15 and tC=10 one can operate the HP. With tH=10°C not if 5°C is lowest level
                - sum(HP_m3Ps[tC,tl+temp_levels.step,t]*3600/(1/verl-1/cop(tC,tl+temp_levels.step,0.5)) for tC in temp_levels if temp_levels.start<tC<tl+temp_levels.step)
            else 0 end
            +
            #Wenn tH=100 ist dann wird temp_levels.stop (95°C) auf die VL + 5°C aufgeheizt und dann in einem Wärmetauscher auf RL+5°C abgekühlt und geht wieder in den SPeicher  
            # m_ab=m3pS_zu*deltaT_zu/deltaT_ab/(1/verl-1/COP)
            if HP_tH_upperbound[t]==supply_levels.stop && (tl+temp_levels.step) == supply_levels.stop
                - sum(HP_m3Ps[tC,supply_levels.stop,t]*3600*temp_levels.step/(VL[t]+temp_levels.step-temp_levels.stop)/(1/verl-1/cop(tC,VL[t]+temp_levels.step,0.5)) for tC in temp_levels if temp_levels.start<tC<supply_levels.stop)
            else 0 end
            )
            )
    # Q_WP_ab_MW=m3pS_zu*deltaT_zu*4.180/(1/verl-1/COP)
    @expression(m,heatPumpsSupplyatSupplytemp_MW[ t in time_steps],
    if VL[t]+temp_levels.step> temp_levels.stop
        sum(HP_m3Ps[tC,supply_levels.stop,t]*temp_levels.step*4.180/(1/verl-1/cop(tC,VL[t]+temp_levels.step,0.5)) for tC in temp_levels if temp_levels.start<tC) 
    else 0 end
    )


    # Pelec_MW=HP_m3Ps[tC,tH,t]*deltaT_zu*4.180/(COP/verl-1)   
    @expression(m,powerCons_HP[ t in time_steps],
        sum(HP_m3Ps[tC,tH,t]*temp_levels.step*4.180/(cop(tC,tH,0.5)/verl-1) for tC in temp_levels for tH in supply_levels if (temp_levels.start<tC<tH<supply_levels.stop)) 
        + 
        if HP_tH_upperbound[t]==supply_levels.stop
                sum(HP_m3Ps[tC,supply_levels.stop,t]*temp_levels.step*4.180/(cop(tC,VL[t]+temp_levels.step,0.5)/verl-1) for tC in temp_levels if temp_levels.start<tC<supply_levels.stop) 
        else 0 end
    )

    @expression(m,powerCons_HP_inStor[ t in time_steps],
    sum(HP_m3Ps[tC,tH,t]*temp_levels.step*4.180/(cop(tC,tH,0.5)/verl-1) for tC in temp_levels for tH in supply_levels if (temp_levels.start<tC<tH<supply_levels.stop)) 
    )

  

    @expression(m,powerCons_HP_forVL[ t in time_steps],
    if HP_tH_upperbound[t]==supply_levels.stop
            sum(HP_m3Ps[tC,supply_levels.stop,t]*temp_levels.step*4.180/(cop(tC,VL[t]+temp_levels.step,0.5)/verl-1) for tC in temp_levels if temp_levels.start<tC<supply_levels.stop) 
    else 0 end
    )


    # Pelec_MW=HP_m3Ps[tC,tH,t]*deltaT_zu*4.180/(COP/verl-1)   
    @constraint(m,hp_Cap_Lim[t in time_steps],powerCons_HP[t]<=HP_Cap_MW )


    @expression(m,totalCost_HP_exp,
    +HP_Cap_MW*annualizedCost_EuropMWpYr(invCostHP_EuroPkW_elec,depreciationYearsHP,wacc)
    )
    # We assume in all HP Modes the same electrical input capacity. If not one would have sum(Hp_Elec_MW(modes)/cap_elec_MW[modes]<=1 for modes in Modes)
    return m
    #=
    storageLevBal_HP_exp          =m[:storageLevBal_HP_exp]
    heatPumpsSupplyatSupplytemp_MW=m[:heatPumpsSupplyatSupplytemp_MW]
    powerCons_HP                  =m[:powerCons_HP]
    totalCost_HP_exp              =m[:totalCost_HP_exp]        
    =#
end

"""
    add_pumpedHydroInCavern_var!(m::Model, temp_levels::OrdinalRange{Int64, Int64}, time_steps::OrdinalRange{Int64, Int64})

DOCSTRING

# Arguments:
- `m`: DESCRIPTION
- `temp_levels`: DESCRIPTION
- `time_steps`: DESCRIPTION
"""
function add_pumpedHydroInCavern_var!(m::Model,temp_levels::OrdinalRange{Int64,Int64},time_steps::OrdinalRange{Int64,Int64}) 
    @variable(m, turbine_MW[temp_levels,time_steps]>=0)
    @variable(m, pumping_MW[temp_levels,time_steps]>=0)
    @variable(m, 10 <=turbineCap_MW)

    #=  
    turbine_MW    = m[:turbine_MW]
    pumping_MW    = m[:pumping_MW]
    turbineCap_MW = m[:turbineCap_MW] 
    =#
    return m
end


function add_pumpedHydroInCavern_exp!(m::Model,cav_levels::Tuple{String,String}, fallofheight_m::Int64, turb_eff::Float64, 
    invCostTurbine_EuroPkW_elec::Float64, depreciationYearsTurbine::Int64, wacc::Float64, temp_levels::OrdinalRange{Int64,Int64},time_steps::OrdinalRange{Int64,Int64})  
    turbine_MW    = m[:turbine_MW]
    pumping_MW    = m[:pumping_MW]
    turbineCap_MW = m[:turbineCap_MW] 
    
    @expression(m,powerProdHyp_exp[ t in time_steps],
        sum(turbine_MW[tl,t] -pumping_MW[tl,t] for tl in temp_levels) 
    )

    @expression(m,storageLevBalHyp_exp[ca in cav_levels::Tuple{String,String},tl in temp_levels::OrdinalRange{Int64,Int64}, t in time_steps],
        (
        if ca=="hi"
            -(turbine_MW[tl,t]*1000/9.81/fallofheight_m/turb_eff -pumping_MW[tl,t]*1000/9.81/fallofheight_m*turb_eff )*3600
        else
            (turbine_MW[tl,t]*1000/9.81/fallofheight_m/turb_eff -pumping_MW[tl,t]*1000/9.81/fallofheight_m*turb_eff )*3600    
        end)

    )
    
    @expression(m,totalCostHyd_exp, turbineCap_MW*annualizedCost_EuropMWpYr(invCostTurbine_EuroPkW_elec,depreciationYearsTurbine,wacc))
    
    #= 
    powerProdHyp_exp     =m[:powerProdHyp_exp] 
    storageLevBalHyp_exp =m[:storageLevBalHyp_exp]
    totalCostHyd_exp     =m[:totalCostHyd_exp]
    =#
    return m
end

function add_pumpedHydroInCavern_con!(m::Model,cav_levels::Tuple{String,String}, fallofheight_m::Int64, turb_eff::Float64,temp_levels::OrdinalRange{Int64,Int64},time_steps::OrdinalRange{Int64,Int64}) 
    turbine_MW    = m[:turbine_MW]
    pumping_MW    = m[:pumping_MW]
    turbineCap_MW = m[:turbineCap_MW] 
    
    @constraint(m, turbineCapLimit_MW[t in time_steps],
        sum(turbine_MW[tl,t] + pumping_MW[tl,t] for tl in temp_levels)<=turbineCap_MW 
    )

    return m
end

function add_wasteHeat_varCon!(m::Model,AbwaermePotential_MW::Array{Float64,1}, temp_levels::OrdinalRange{Int64,Int64},time_steps::OrdinalRange{Int64,Int64}) 

    heatExchangerCap_Wasteheat_MW=m[:heatExchangerCap_Wasteheat_MW]

    # Waste heat wird mit 95°C geliefert nach Wärmetauscher nur mehr 90°C
    @variable(m, 0<=wasteheat_90C_m3[t in time_steps])#<=AbwaermePotential_MW[t]

    # return flow of cool water to wasteheat factory
    @variable(m, returnWasteH_m3[tl in temp_levels,time_steps;tl <90]>=0)

    @constraint(m, returnWasteHBal[t in time_steps],sum(returnWasteH_m3[tl,t] for tl in temp_levels if tl<90)==wasteheat_90C_m3[t])

    @constraint(m, wasteHeatPotLim[t in time_steps], sum( returnWasteH_m3[tl,t]*(90-tl)*4.180 for tl in temp_levels if tl<90) <=AbwaermePotential_MW[t])

    @constraint(m, heatExCap[t in time_steps], sum( returnWasteH_m3[tl,t]*(90-tl)*4.180 for tl in temp_levels if tl<90) <=heatExchangerCap_Wasteheat_MW)

    #= 
    wasteheat_90C_m3=m[:wasteheat_90C_m3]
    returnWasteH_m3=m[:returnWasteH_m3]
    =#
    return m
end


function buildModel(m::Model,cav_levels::Tuple{String,String},temp_levels::OrdinalRange{Int64,Int64},supply_levels::OrdinalRange{Int64,Int64},
    powerprice::Array{Float64,1},districtHeatPrice::Array{Float64,1},AbwaermePotential_MW::Array{Float64,1},
    VL::Array{Float64,1}, HP_tH_upperbound::Array{Float64,1}, RL::Array{Float64,1},
    invCostCavern_EuroPm3::Float64        ,depreciationYearsCavern::Int64,
    invCostHP_EuroPkW_elec::Float64       ,depreciationYearsHP::Int64,
    invCostTurbine_EuroPkW_elec::Float64  ,depreciationYearsTurbine::Int64,
    costHeatExchanger_wasteheat_EuroPkWth::Float64,depreciationYearsHeatEx::Int64,
    CapLimitFW_Waermetauscher_MW::Float64,
    CapLimitWasteH_Waermetauscher_MW::Float64,
    CapLimitRiver_Waermetauscher_MW::Float64,
    riverTemp::Array{Float64,1},
    wacc::Float64,
    capOneCavern_1000m3::Int64, fallofheight_m::Int64, turb_eff::Float64,
    time_steps::OrdinalRange{Int64,Int64})


    add_cavern_var!(m,cav_levels,temp_levels,VL, CapLimitFW_Waermetauscher_MW,CapLimitWasteH_Waermetauscher_MW,CapLimitRiver_Waermetauscher_MW,riverTemp,time_steps)
    add_heatpumpInCavern_var!(m,temp_levels,supply_levels,HP_tH_upperbound,time_steps)
    add_pumpedHydroInCavern_var!(m,temp_levels,time_steps)
    add_wasteHeat_varCon!(m,AbwaermePotential_MW,temp_levels,time_steps) 
    
    add_heatpumpInCavern_expCon!(m,temp_levels,supply_levels, time_steps,VL, HP_tH_upperbound, RL, invCostHP_EuroPkW_elec,depreciationYearsHP, wacc)
    add_pumpedHydroInCavern_exp!(m,cav_levels, fallofheight_m, turb_eff, invCostTurbine_EuroPkW_elec,depreciationYearsTurbine, wacc, temp_levels,time_steps)   
    
    
    add_cavern_conUexp!(m::Model,cav_levels,capOneCavern_1000m3,invCostCavern_EuroPm3,depreciationYearsCavern,costHeatExchanger_wasteheat_EuroPkWth,depreciationYearsHeatEx, wacc, temp_levels,time_steps,binary_period)
    add_pumpedHydroInCavern_con!(m,cav_levels, fallofheight_m, turb_eff,temp_levels,time_steps) 
 
    storageLevBal_HP_exp          =m[:storageLevBal_HP_exp]
    heatPumpsSupplyatSupplytemp_MW=m[:heatPumpsSupplyatSupplytemp_MW]
    powerCons_HP                  =m[:powerCons_HP]
    totalCost_HP_exp              =m[:totalCost_HP_exp]        
    powerProdHyp_exp     =m[:powerProdHyp_exp] 
    storageLevBalHyp_exp =m[:storageLevBalHyp_exp]
    totalCostHyd_exp     =m[:totalCostHyd_exp]
    M3_storage=m[:M3_storage]
    totalCostCavern_exp=m[:totalCostCavern_exp]
    wasteheat_90C_m3=m[:wasteheat_90C_m3]
    returnWasteH_m3=m[:returnWasteH_m3]
    districtHeatFromStorage_m3pS=m[:districtHeatFromStorage_m3pS]
    heatExchangerCap_CavToFW_MW=m[:heatExchangerCap_CavToFW_MW]
    returnRiverwater_m3pS=m[:returnRiverwater_m3pS]
    riverwater_m3pS=m[:riverwater_m3pS]
    heatExchangerCap_river_MW=m[:heatExchangerCap_river_MW]
    
    @expression(m, elecCons_Exp[t in time_steps],powerCons_HP[t]-powerProdHyp_exp[t])

    #expression Mone from pumping, heapump, districtheat fro stprage district heat for suply at supl


    # design point is 95-->65,60-->90 mit t_m=5
    @constraint(m,heatExchangerCap_Cavern_Limit[t in time_steps],
        # the district heat can be also supplied from templevels above VL+temp_levels.step, 
        # but then the t_m is increased and therefore also the power of the heatexchanger
        # The mass flow is always adjusted such that the flow is cooled to RL[t]+temp_levels.step
        (if (VL[t]+temp_levels.step<=temp_levels.stop) 
            sum(districtHeatFromStorage_m3pS[tl,t]*(tl-RL[t]-temp_levels.step)*4.180/t_m(tl-VL[t],5)*5 for tl in temp_levels if tl>=VL[t]+temp_levels.step)
        else 0 end   
        +heatPumpsSupplyatSupplytemp_MW[t]
        )<=heatExchangerCap_CavToFW_MW)


    @constraint(m,heatExchangerCap_river_Limit[t in time_steps; riverTemp[t]>=(temp_levels.start+2*temp_levels.step)],
                sum(returnRiverwater_m3pS[tl,t]*(riverTemp[t]-temp_levels.step-tl)*4.180 for tl in temp_levels if riverTemp[t]>=(tl+2*temp_levels.step))
                <=heatExchangerCap_river_MW)

        #@constraint(m,riverBal[t in time_steps;
        #riverTemp[t]>=(temp_levels.start+2*temp_levels.step)],sum(returnRiverwater_m3pS[tl,t] for tl in temp_levels if (tl+2*temp_levels.step)<=riverTemp[t])==riverwater_m3pS[t])
    

    # Q_MW=m3pS*deltaT*4.180
    @constraint(m, storageLevelBalance_con[h in cav_levels,tl in temp_levels, t in time_steps], 
        (if(h=="hi")  
            (M3_storage[h, tl, t] +storageLevBalHyp_exp[h,tl,t] +storageLevBal_HP_exp[tl,t] 
            + if tl==90; wasteheat_90C_m3[t]*3600 else 0 end 
            - if tl<90; returnWasteH_m3[tl,t]*3600  else 0 end 
            + if (tl+temp_levels.step)==riverTemp[t]>=(temp_levels.start+2*temp_levels.step); riverwater_m3pS[t]*3600 else 0 end
            - if (tl+2*temp_levels.step)<=riverTemp[t]; returnRiverwater_m3pS[tl,t]*3600 else 0 end
            
            - if tl>=VL[t]+temp_levels.step;  districtHeatFromStorage_m3pS[tl,t]*3600 else 0 end 
            # districtHeatFromStorage_m3pS[t] exists only if VL[t]+temp_levels.step<=temp_levels.stop
            + if ((tl==RL[t]+temp_levels.step) && (VL[t]+temp_levels.step<=temp_levels.stop));  sum(districtHeatFromStorage_m3pS[tl2,t] for tl2 in temp_levels if tl2>=(VL[t]+temp_levels.step))*3600 else 0 end 
            )
        else
            M3_storage[h, tl, t]+ storageLevBalHyp_exp[h,tl,t]
        end)
        ==M3_storage[h, tl, rem(t,time_steps.stop)+1]
    )

    @objective(m, Min, 
        (totalCostCavern_exp +totalCost_HP_exp +totalCostHyd_exp +
        sum((powerCons_HP[t]-powerProdHyp_exp[t])*powerprice[t] 
        - districtHeatPrice[t]*
            (heatPumpsSupplyatSupplytemp_MW[t]
            + sum(districtHeatFromStorage_m3pS[tl,t]*(tl-RL[t]-temp_levels.step)*4.180 for tl in temp_levels if tl>=(VL[t]+temp_levels.step))
            )
        for t in time_steps)
        )
    )

end


time_steps=(1:8760)
binary_period=(1:0)  # Keine Binaries: binary_period=(1:0)
temp_levels=(5:5:95)
#Info:  HP_m3Ps[tC,tH,t] ist the m3/s mass flow on the source side and takes source heat from tC and cools it to tC-5 and supplys heat at tH by reheating tH-5
#With tH=15 and tC=10 one can operate the HP. With tH=10°C not if 5°C is lowest level
supply_levels=(15:5:100)# 100 is a flag value to indicate the supply temperature is the one of the district heating grid 


A=XLSX.readdata("Input.xlsx", "Tabelle1", "B5:J8765")
districtHeatPrice=broadcast(max,convert(Array{Float64,1},A[2:8761,2]),0) # Euro/MWh
powerprice=convert(Array{Float64,1},A[2:8761,3]) # Euro/MWh
Aussentemp_C=convert(Array{Float64,1},A[2:8761,4])
VL=convert(Array{Float64,1},A[2:8761,5])
RL=convert(Array{Float64,1},A[2:8761,6])
FW_Demand_MW=convert(Array{Float64,1},A[2:8761,7])
AbwaermePotential_MW=convert(Array{Float64,1},A[2:8761,8])
riverTemp=round.(convert(Array{Float64,1},A[2:8761,9])./5).*5

cav_levels=("hi","lo")
capOneCavern_1000m3=10
fallofheight_m=700
turb_eff=0.9
invCostCavern_EuroPm3=140.0 #140 Euro/m3
invCostHP_EuroPkW_elec=200.0   #200 Euro/kW_elec: 600 Euro/kWth mit Cop =3 ergibt 200 Euro/kW_elec
#C:\Gerhard\Model\Julia\SeasonalStorage\HiREPS Originals\Kosten und Umrechnung Kronen in Euro.xlsx
invCostTurbine_EuroPkW_elec=1000.0   #1000 Euro/kW
depreciationYearsCavern=60
depreciationYearsHP=15
depreciationYearsTurbine=35
depreciationYearsHeatEx=20
CapLimitFW_Waermetauscher_MW = 70.0 # MW
CapLimitWasteH_Waermetauscher_MW=10.5 # MW
CapLimitRiver_Waermetauscher_MW = 200.0 # MW
#costHeatExchanger_groundwater=272# Euro/kWthermal_output
#costHeatExchanger_fluegasCondensation=101# Euro/kWthermal_output
costHeatExchanger_wasteheat_EuroPkWth=85.0 # Euro/kWthermal_output 

verl=0.97
wacc=5.0

# round VL und RL to 5 degree steps
for i in eachindex(VL)
    VL[i]=
        # Bis 100 in 5° Steps runden darüber kontinuierlich
        if VL[i] > temp_levels.stop+temp_levels.step
            VL[i]
        else
            round(VL[i]/5)*5
        end
end

for i in eachindex(RL)
    RL[i]=round(RL[i]/5)*5
end

HP_tH_upperbound=convert(Array{Float64,1},[(VL[t]+temp_levels.step>temp_levels.stop) ? supply_levels.stop : temp_levels.stop for t in time_steps ])


buildModel(m,cav_levels,temp_levels,supply_levels,powerprice,districtHeatPrice,AbwaermePotential_MW,VL, HP_tH_upperbound, RL,
    invCostCavern_EuroPm3,depreciationYearsCavern,
    invCostHP_EuroPkW_elec,depreciationYearsHP,
    invCostTurbine_EuroPkW_elec, depreciationYearsTurbine,
    costHeatExchanger_wasteheat_EuroPkWth,depreciationYearsHeatEx,
    CapLimitFW_Waermetauscher_MW,CapLimitWasteH_Waermetauscher_MW,CapLimitRiver_Waermetauscher_MW,
    riverTemp,
    wacc, capOneCavern_1000m3, fallofheight_m, turb_eff, time_steps)

    
JuMP.optimize!(m)

termination_status(m) # why did the solver stop https://jump.dev/JuMP.jl/stable/solutions/#Termination-statuses-1

primal_status(m) #https://jump.dev/JuMP.jl/stable/solutions/#Solution-statuses-1

dual_status(m) #https://jump.dev/JuMP.jl/stable/solutions/#Solution-statuses-1

objective_value(m)


#JuMP.write_to_file(m, "SeasonalM.lp")
#CPXwriteprob( mm.env, mm.lp, "mymodel2.lp", "LP" )
#stat = CPXgetstat(mm.env, mm.lp)
#print(stat)


if termination_status(m) == MOI.OPTIMAL
	@show objective_value(m);
	print("") # damit der Wert nicht 2 mal geprintet wird
end
