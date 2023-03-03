/** Abstraction for a Reaction */
class Reaction {
    constructor( libsbmlReaction, functionDefinitions, libsbml, pyodide ) {
        this.reaction = libsbmlReaction;
        this.reactants = [];
        this.products = [];
        var i;
        for (i = 0; i < this.reaction.getNumReactants(); i++) {
            this.reactants.push(this.reaction.getReactant(i));
        }
        for (i = 0; i < this.reaction.getNumProducts(); i++) {
            this.products.push(this.reaction.getProduct(i));
        }
        this.id = this.reaction.getId();
        this.kineticLaw = new KineticLaw(this.reaction.getKineticLaw(), this, functionDefinitions, libsbml, pyodide);
    }

    getId() {
        return this.id;
    }

    toString() {
        var i;
        var reactantStr = "";
        var productStr = "";
        for (i = 0; i < this.reactants.length - 1; i++) {
            reactantStr += `${this.reactants[i].getSpecies()} + `;
        }
        for (i = 0; i < this.products.length - 1; i++) {
            productStr += `${this.products[i].getSpecies()} + `;
        }
        if (this.reactants.length > 0) {
            reactantStr += this.reactants[this.reactants.length - 1].getSpecies();
        }
        if (this.products.length > 0) {
            productStr += this.products[this.products.length - 1].getSpecies();
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