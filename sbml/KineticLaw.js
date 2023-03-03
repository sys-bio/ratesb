/** Provides Information on SBML Kinetics Laws */

class KineticLaw {
    constructor( libsbmlKinetics, reaction, functionDefinitions, libsbml, pyodide) {
        this.libsbmlKinetics = libsbmlKinetics;
        this.formula = this.libsbmlKinetics.getFormula();
        this.reaction = reaction;
        this._curDepth = 0;
        try {
            this.symbols = this._getSymbols();
        } catch (e) {
            this.symbols = [];
        }
        if (functionDefinitions == None) {
            this.expandedFormula = null;
        } else {
            this.expandFormula(functionDefinitions);
        }
        this.expressionFormula = null;
    }

    expandFormula(functionDefinitions) {
        var maxRecursion = 5;
        this.expandedFormula = this._expandFormula(this.formula, functionDefinitions, 0, maxRecursion);
    }

    _expandFormula(expansion, functionDefinitions, numRecursion, maxRecursion) {
        if (numRecursion > maxRecursion) {
            return expansion;
        }
        var done = true;
        for (var i = 0; i < functionDefinitions.length; i++) {
            var fd = functionDefinitions[i];
            var rx = new RegExp(String.raw(`${fd.id}\(.*?\)`));
            var matches = [];
            var match;
            while ((match = rx.exec(expansion)) !== null) {
                matches.push(match);
            }
            if (matches.length == 0) {
                continue
            }
            done = false;
            for (var j = 0; j < matches.length; j++) {
                
            }
        }
    }

    
    _augment(astNode, result, maxDepth) {
        this._curDepth++;
        if (this._curDepth > maxDepth) {
            throw new Error("Bad Kinetics Math");
        }
        for (var i = 0; i < astNode.getNumChildren(); i++) {
            var childNode = astNode.getChild(i);
            var additions;
            if (childNode.getName() == null) {
                additions = this._augment(childNode, result, maxDepth);
                result.push(additions);
            } else {
                if (childNode.isFunction()) {
                    additions = this._augment(childNode, result, maxDepth);
                    result.push(additions);
                } else {
                    result.push(childNode.getName());
                }
            }
        }
        return result;
    }

    _getSymbols() {
        var maxDepth = 20;
        this._curDepth = 0;
        var astNode = this.libsbmlKinetics.getMath();
        var result;
        if (astNode.getName() == null) {
            result = [];
        } else {
            result = [astNode.getName()];
        }
        return this._augment(astNode, result, maxDepth);
    }

    

}