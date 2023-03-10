/** Process SBML file with an input SBML model in string */
class ProcessSBML {
    constructor( modelStr, pyodide, libsbml ) {
        this.pyodide = pyodide;
        this.libsbml = libsbml;
        const reader = new this.libsbml.SBMLReader();
        const sbmlDoc = reader.readSBMLFromString(modelStr);
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
        this.functionDefinitions = this.getFunctionDefinitions();
        this.reactions = [];
        for (i = 0; i < this.model.getNumReactions(); i++) {
            this.reactions.push(new Reaction(this.model.getReaction(i), this.functionDefinitions, libsbml, pyodide));
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
        const species = new Array(this.model.getNumSpecies());
        for (var i = 0; i < this.model.getNumSpecies(); i++) {
          species[i] = this.model.getSpecies(i).getId();
        }
        return species;
    }

    getParameterList() {
        const parameters = new Array(this.model.getNumParameters());
        for (var i = 0; i < this.model.getNumParameters(); i++) {
          parameters[i] = this.model.getParameter(i).getId();
        }
        return parameters;
    }
}