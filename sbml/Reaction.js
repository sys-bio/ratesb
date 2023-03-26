/** Abstraction for a Reaction */
class Reaction {
    constructor( libsbmlModel, libsbmlReaction, functionDefinitions, pyodide, processSBML ) {
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
        this.kineticLaw = new KineticLaw(libsbmlModel, libsbmlReaction, this.reaction.getKineticLaw(), this, functionDefinitions, pyodide, processSBML);
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
        return `${reactantStr} -> ${productStr}; ${kineticStr}`;
    }
}