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
        if (!(this.kineticLaw.classificationCp["zerothOrder"]
            || this.kineticLaw.classificationCp["powerTerms"]
            || this.kineticLaw.classificationCp["UNDR"]
            || this.kineticLaw.classificationCp["UNMO"]
            || this.kineticLaw.classificationCp["BIDR"]
            || this.kineticLaw.classificationCp["BIMO"]
            || this.kineticLaw.classificationCp["MM"]
            || this.kineticLaw.classificationCp["MMcat"]
            || this.kineticLaw.classificationCp["Hill"]
            || this.kineticLaw.classificationCp["Polynomial"]
            )) {
            ret += `\nWarning: Unrecognized rate law from the standard list.`;
        }
        if (this.kineticLaw.analysis["floatingSpecies"].length > 0) {
            ret += `\nWarning: floating reactant: ${this.kineticLaw.analysis["floatingSpecies"].toString()}`;
        }
        if (this.kineticLaw.analysis["inconsistentProducts"].length > 0) {
            ret += `\nWarning: irreversible reaction kinetic law contains products: ${this.kineticLaw.analysis["inconsistentProducts"].toString()}`;
        }
        if (this.kineticLaw.analysis["nonIncreasingSpecies"].length > 0) {
            ret += `\nWarning: non increasing species: ${this.kineticLaw.analysis["nonIncreasingSpecies"].toString()}`;
        }
        if (this.kineticLaw.analysis["namingConvention"]['k'].length > 0) {
            ret += `\nWarning: elements should start with k: ${this.kineticLaw.analysis["namingConvention"]['k'].toString()}`;
        }
        if (this.kineticLaw.analysis["namingConvention"]['v'].length > 0) {
            ret += `\nWarning: elements should start with v: ${this.kineticLaw.analysis["namingConvention"]['v'].toString()}`;
        }
        return ret;
    }
}