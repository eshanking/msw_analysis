package DoseResponse;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
// import HAL.lib.CommandLine.Help.TextTable.Cell;


class DoseResponseCell extends AgentSQ2Dunstackable<DoseResponseGrid> {
    int genotype;
    double g_drugless;
    double ic50;
    double hillCoef;
    double gmin;

    // public DoseResponseCell(int genotype){
    //     this.genotype = genotype;
    //     this.g_drugless = G.GrowthRateList[genotype];
    //     this.ic50 = G.ic50List[genotype];
    //     this.hillCoef = G.hillCoefList[genotype];
    //     this.gmin = G.gminList[genotype];
    // }

    public int GetGenotype(){
        return this.genotype;
    }



    public double HillEqn(double conc){
        // convert conc to log10
        // conc = conc*G.drugConcScale;
        
        if (G.threeParamHill){ // 3 parameter hill equation
            return this.g_drugless / (1 + Math.exp((this.ic50 - Math.log10(conc)) / this.hillCoef));
        }
        else{ // 4 parameter hill equation
            // if (conc < Math.pow(10,-5)){ // prevent underflow error
            //     return this.g_drugless;
            // }
            // else{
            //     return this.g_drugless + ((this.gmin - this.g_drugless) * Math.pow(conc, this.hillCoef)) / (Math.pow(this.ic50, this.hillCoef) + Math.pow(conc, this.hillCoef));
            // }
            return this.g_drugless + ((this.gmin - this.g_drugless) * Math.pow(conc, this.hillCoef)) / (Math.pow(this.ic50, this.hillCoef) + Math.pow(conc, this.hillCoef));
        }
    }

    public double GetFitness(double conc){
        double scaledConc = conc*G.drugConcScale;
        if (G.useMaxConc){
            if (scaledConc > G.maxConc){
                return 0;
            } else{
                return HillEqn(scaledConc);
            }
        }
        else{
            return HillEqn(scaledConc);
        }
    }

    public void setParams(int genotype){
        this.genotype = genotype;
        this.g_drugless = G.GrowthRateList[genotype];
        this.ic50 = G.ic50List[genotype];
        this.hillCoef = G.hillCoefList[genotype];
        this.gmin = G.gminList[genotype];
    
    }
}

    // public double hillFn(double conc, double gmax, double gmin, double hc, double ic50) {
    //     double y = gmax + ((gmin - gmax) * Math.pow(conc, hc)) / (Math.pow(ic50, hc) + Math.pow(conc, hc));
    //     return y;
    // }