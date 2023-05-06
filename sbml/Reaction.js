/** Abstraction for a Reaction */
class Reaction {
    constructor( libsbmlModel, libsbmlReaction, functionDefinitions, pyodide, processSBML, namingConvention ) {
        this.reaction = libsbmlReaction;
        console.log(this.reaction)
        this.reactantList = [];
        var i;
        for (i = 0; i < this.reaction.getNumReactants(); i++) {
            var reactant = this.reaction.getReactant(i);
            this.reactantList.push(reactant.getSpecies());
        }
        this.productList = [];
        for (i = 0; i < this.reaction.getNumProducts(); i++) {
            var product = this.reaction.getProduct(i)
            this.productList.push(product.getSpecies());
        }
        this.id = this.reaction.getId();
        this.namingConvention = namingConvention;
        this.kineticLaw = new KineticLaw(libsbmlModel, libsbmlReaction, this.reaction.getKineticLaw(), this, functionDefinitions, pyodide, processSBML, namingConvention);
    }

    getId() {
        return this.id;
    }

    toString() {
        var i;
        var reactantStr = "";
        var productStr = "";
        for (i = 0; i < this.reactantList.length - 1; i++) {
            reactantStr += `${this.reactantList[i]} + `;
        }
        for (i = 0; i < this.productList.length - 1; i++) {
            productStr += `${this.productList[i]} + `;
        }
        if (this.reactantList.length > 0) {
            reactantStr += this.reactantList[this.reactantList.length - 1];
        }
        if (this.productList.length > 0) {
            productStr += this.productList[this.productList.length - 1];
        }
        var kineticStr = "";
        if (this.kineticLaw.expandedFormula !== null) {
            kineticStr = this.kineticLaw.expandedFormula;
        } else {
            kineticStr = this.kineticLaw.formula;
        }
        var ret = `${reactantStr} -> ${productStr}; ${kineticStr}`;
        if (kineticLaw.analysis["emptyRateLaw"]) {
            ret += "\nError E001: No Rate Law entered";
        } else if (kineticLaw.analysis["pureNumber"]) {
            ret += "\nWarning W001: Rate law contains only number";
        } else {
            if (kineticLaw.analysis["floatingSpecies"].length > 0) {
                ret += `\nError E002: Floating reactant: ${kineticLaw.analysis["floatingSpecies"].toString()}`;
            }
            if (!(
                kineticLaw.classificationCp["zerothOrder"] ||
                kineticLaw.classificationCp["powerTerms"] ||
                kineticLaw.classificationCp["UNDR"] ||
                kineticLaw.classificationCp["UNMO"] ||
                kineticLaw.classificationCp["BIDR"] ||
                kineticLaw.classificationCp["BIMO"] ||
                kineticLaw.classificationCp["MM"] ||
                kineticLaw.classificationCp["MMcat"] ||
                kineticLaw.classificationCp["Hill"] ||
                kineticLaw.classificationCp["Polynomial"]
            )) {
                ret += "\nWarning W002: Unrecognized rate law from the standard list.";
            }
            if (kineticLaw.analysis["inconsistentProducts"].length > 0) {
                ret += `\nWarning W004: Irreversible reaction kinetic law contains products: ${kineticLaw.analysis["inconsistentProducts"].toString()}`;
            }
            if (kineticLaw.analysis["nonIncreasingSpecies"].length > 0) {
                ret += `\nWarning W005: Non increasing species: ${kineticLaw.analysis["nonIncreasingSpecies"].toString()}`;
            }
            if (kineticLaw.analysis["namingConvention"]['k'].length > 0) {
                ret += `\nWarning W006: We recommend that these parameters start with 'k': ${kineticLaw.analysis["namingConvention"]['k'].toString()}`;
            }
            if (kineticLaw.analysis["namingConvention"]['v'].length > 0) {
                ret += `\nWarning W007: We recommend that these parameters start with 'v': ${kineticLaw.analysis["namingConvention"]['v'].toString()}`;
            }
            if (kineticLaw.analysis["formattingConvention"] == 1) {
                ret += "\nWarning W008: Formatting convention not followed (compartment before parameters before species).";
            }
            if (kineticLaw.analysis["formattingConvention"] == 2) {
                ret += "\nWarning W009: Elements of the same type should be ordered alphabetically.";
            }
        }
        return ret;
    }
}