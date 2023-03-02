/** Abstraction for a Reaction */
class Reaction {
    constructor( libsbmlReaction, functionDefinitions ) {
        this.reaction = libsbmlReaction;
        this.reactants = [];
    }
}

class FunctionDefinition {
    constructor( sbmlFunctionDefinition, libsbml) {
        this.libsbml = libsbml;
        this.sbmlFunctionDefinition = sbmlFunctionDefinition;
        this.functionName = this.sbmlFunctionDefinition.getName();
        this.id = this.sbmlFunctionDefinition.getId();
        this.argumentNames = []
        for (var i = 0; i < this.sbmlFunctionDefinition.getNumArguments(); i++) {
            this.argumentNames.push(this.sbmlFunctionDefinition.getArgument(i).getName());
        }
        this.body = new libsbml.SBMLFormulaParser().formulaToL3String(this.sbmlFunctionDefinition.getBody());
    }

    toString() {
        argumentCall = this.argumentNames.toString();
        callStr = `${this.id}(${argumentCall})`;
        return `${callStr}: ${this.body}`;
    }
}

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
        console.log(this.functionDefinitions);
        // this.reactions = [];
        // for (i = 0; i < this.model.getNumReactions(); i++) {
        //     this.species.push(this.model.getReaction(i), this.functionDefinitions);
        // }
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