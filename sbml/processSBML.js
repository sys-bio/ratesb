/** Process SBML file with an input SBML model in string */
class ProcessSBML {
    constructor( modelStr, pyodide, libsbml, namingConvention ) {
        this.pyodide = pyodide;
        this.libsbml = libsbml;
        const reader = new this.libsbml.SBMLReader();
        const sbmlDoc = reader.readSBMLFromString(modelStr);
        // const converter = libsbml.LocalParameter;
        // console.log(converter)
        // converter.setDocument(sbmlDoc);
        // converter.convert();
        const libsbmlModel = sbmlDoc.getModel();
        this.model = libsbmlModel;
        console.log(this.model);
        this.reactions = [];
        this.species = [];
        this.parameters = [];
        this.functionDefinitions = [];
        var i;
        for (i = 0; i < this.model.getNumSpecies(); i++) {
            this.species.push(this.model.getSpecies(i));
        }
        for (i = 0; i < this.model.getNumParameters(); i++) {
            this.parameters.push(this.model.getParameter(i));
        }
        this.speciesList = this.getSpeciesList();
        this.parameterList = this.getParameterList();
        this.compartmentList = this.getCompartmentList();
        console.log(this.speciesList)
        if (this.speciesList.length > 0) {
            this.pyodide.runPython(`
                ${this.speciesList.toString()} = sympy.symbols("${this.speciesList.join(" ")}")
            `);
        }
        console.log(this.parameterList)
        if (this.parameterList.length > 0) {
            this.pyodide.runPython(`
                ${this.parameterList.toString()} = sympy.symbols("${this.parameterList.join(" ")}")
            `);
        }
        console.log(this.compartmentList)
        if (this.compartmentList.length > 0) {
            this.pyodide.runPython(`
                ${this.compartmentList.toString()} = sympy.symbols("${this.compartmentList.join(" ")}")
            `);
        }
        this.functionDefinitions = this.getFunctionDefinitions();
        this.reactions = [];
        for (i = 0; i < this.model.getNumReactions(); i++) {
            this.reactions.push(new Reaction(this.model, this.model.getReaction(i), this.functionDefinitions, pyodide, this, namingConvention));
        }
    }

    getFunctionDefinitions() {
        var functionDefinitions = [];
        for (var i = 0; i < this.model.getNumFunctionDefinitions(); i++) {
            functionDefinitions.push(new FunctionDefinition(this.model.getFunctionDefinition(i), this.libsbml));
        }
        return functionDefinitions;
    }

    // Get all species in the model
    getSpeciesList() {
        const species = [];
        for (var i = 0; i < this.model.getNumSpecies(); i++) {
          species.push(this.model.getSpecies(i).getId());
        }
        return species;
    }

    getParameterList() {
        const parameters = [];
        for (var i = 0; i < this.model.getNumParameters(); i++) {
          parameters.push(this.model.getParameter(i).getId());
        }
        return parameters;
    }

    getCompartmentList() {
        const compartments = [];
        for (var i = 0; i < this.model.getNumCompartments(); i++) {
            compartments.push(this.model.getCompartment(i).getId());
        }
        return compartments;
    }
}