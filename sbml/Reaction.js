/** Abstraction for a Reaction */
class Reaction {
    constructor( libsbmlModel, libsbmlReaction, functionDefinitions, pyodide, processSBML, checks ) {
        this.reaction = libsbmlReaction;
        console.log(this.reaction)
        this.boundarySpeciesList = [];

        this.reactantList = [];
        var i;
        for (i = 0; i < this.reaction.getNumReactants(); i++) {
            const reactant = this.reaction.getReactant(i);
            const reactantName = reactant.getSpecies()
            this.reactantList.push(reactantName);
            const libsbmlSpecies = libsbmlModel.getSpecies(reactantName);
            if (libsbmlSpecies.getBoundaryCondition()) {
                this.boundarySpeciesList.push(reactantName);
            }
        }
        this.productList = [];
        for (i = 0; i < this.reaction.getNumProducts(); i++) {
            const product = this.reaction.getProduct(i);
            const productName = product.getSpecies();
            this.productList.push(productName);
            const libsbmlSpecies = libsbmlModel.getSpecies(productName);
            if (libsbmlSpecies.getBoundaryCondition()) {
                this.boundarySpeciesList.push(productName);
            }
        }
        
        this.sortedSpeciesList = this.reactantList.concat(this.productList);
        this.id = this.reaction.getId();
        this.kineticLaw = new KineticLaw(libsbmlModel, libsbmlReaction, this.reaction.getKineticLaw(), this, functionDefinitions, pyodide, processSBML, checks);
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
        var reactionArrow = "=>";
        if (this.reaction.getReversible()) {
            reactionArrow = "->";
        }
        var ret = `${reactantStr} ${reactionArrow} ${productStr}; ${kineticStr}`;
        if (this.kineticLaw.analysis["emptyRateLaw"]) {
            ret += "\nError E001: No Rate Law entered";
        } else if (this.kineticLaw.analysis["pureNumber"]) {
            ret += "\nWarning W001: Rate law contains only number";
        } else {
            if (this.kineticLaw.analysis["floatingSpecies"].length > 0) {
                ret += `\nError E002: Expecting reactant in rate law but not found: ${this.kineticLaw.analysis["floatingSpecies"].toString()}`;
            }
            if (this.kineticLaw.analysis["boundaryFloatingSpecies"].length > 0) {
                ret += `\nWarning W016: Expecting boundary species reactant in rate law but not found: ${this.kineticLaw.analysis["boundaryFloatingSpecies"].toString()}`
            }
            if (!(
                this.kineticLaw.classificationCp["zerothOrder"] ||
                this.kineticLaw.classificationCp["powerTerms"] ||
                this.kineticLaw.classificationCp["UNDR"] ||
                this.kineticLaw.classificationCp["UNMO"] ||
                this.kineticLaw.classificationCp["BIDR"] ||
                this.kineticLaw.classificationCp["BIMO"] ||
                this.kineticLaw.classificationCp["MM"] ||
                this.kineticLaw.classificationCp["MMcat"] ||
                this.kineticLaw.classificationCp["Hill"] ||
                this.kineticLaw.classificationCp["Polynomial"]
            )) {
                ret += "\nWarning W002: Unrecognized rate law from the standard list.";
            }
            if (this.kineticLaw.analysis["inconsistentProducts"].length > 0) {
                ret += `\nWarning W003: Irreversible reaction kinetic law contains products: ${this.kineticLaw.analysis["inconsistentProducts"].toString()}`;
            }
            if (this.kineticLaw.analysis["nonIncreasingSpecies"][0].length > 0) {
                ret += `\nWarning W004: Flux is not increasing as reactant increases: ${this.kineticLaw.analysis["nonIncreasingSpecies"][0].toString()}`;
            }
            if (this.kineticLaw.analysis["nonIncreasingSpecies"][1].length > 0) {
                ret += `\nWarning W0015: Flux is not decreasing as product increases: ${this.kineticLaw.analysis["nonIncreasingSpecies"][1].toString()}`;
            }
            if (this.kineticLaw.analysis["namingConvention"]['k'].length > 0) {
                ret += `\nWarning W005: We recommend that these parameters start with 'k': ${this.kineticLaw.analysis["namingConvention"]['k'].toString()}`;
            }
            if (this.kineticLaw.analysis["namingConvention"]['K'].length > 0) {
                ret += `\nWarning W006: We recommend that these parameters start with 'K': ${this.kineticLaw.analysis["namingConvention"]['k'].toString()}`;
            }
            if (this.kineticLaw.analysis["namingConvention"]['V'].length > 0) {
                ret += `\nWarning W007: We recommend that these parameters start with 'V': ${this.kineticLaw.analysis["namingConvention"]['v'].toString()}`;
            }
            switch (this.kineticLaw.analysis["formattingConvention"]) {
                case 0:
                    // Formatted correctly, no warning needed
                    break;
                case 1:
                    ret += "\nWarning W008: Elements of the same type are not ordered properly.";
                    break;
                case 2:
                    ret += "\nWarning W009: Formatting convention not followed (compartment before parameters before species).";
                    break;
                case 3:
                    ret += "\nWarning W010: Denominator not in alphabetical order.";
                    break;
                case 4:
                    ret += "\nWarning W011: Numerator and denominator not in alphabetical order.";
                    break;
                case 5:
                    ret += "\nWarning W012: Numerator convention not followed and denominator not in alphabetical order.";
                    break;
                case 6:
                    ret += "\nWarning W013: Denominator convention not followed.";
                    break;
                case 7:
                    ret += "\nWarning W014: Numerator not in alphabetical order and denominator convention not followed.";
                    break;
                case 8:
                    ret += "\nWarning W015: Numerator and denominator convention not followed.";
                    break;
                default:
                    // Unhandled case, add an appropriate message if needed
                    break;
            }
        }
        return ret;
    }
}