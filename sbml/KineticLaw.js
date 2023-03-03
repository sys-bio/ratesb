/** Provides Information on SBML Kinetics Laws */

class KineticLaw {
    constructor( libsbmlKinetics, reaction, functionDefinitions, libsbml, pyodide) {
        this.pyodide = pyodide;
        this.libsbmlKinetics = libsbmlKinetics;
        this.formula = this.libsbmlKinetics.getFormula();
        this.reaction = reaction;
        this._curDepth = 0;
        try {
            this.symbols = this._getSymbols();
        } catch (e) {
            this.symbols = [];
        }
        if (functionDefinitions == null) {
            this.expandedFormula = null;
        } else {
            this.expandFormula(functionDefinitions);
        }
        this.expressionFormula = null;
    }

    toString() {
        if (this.expandedFormula == null) {
            return self.formula;
        }
        return self.expandedFormula;
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
            var rx = new RegExp(String.raw`${fd.id}\(.*?\)`);
            var matches = rx.exec(expansion);
            var body = String(fd.body);
            if (matches == null) {
                continue
            }
            done = false;
            for (var j = 0; j < matches.length; j++) {
                var rx_arguments = new RegExp(String.raw`\(.*?\)`);
                var argument_match = rx_arguments.exec(matches[j])[0];
                argument_match = argument_match.trim();
                argument_match = argument_match.substring(1, argument_match.length - 1);
                var call_arguments = argument_match.split(",");
                for (var k = 0; k < call_arguments.length; k++) {
                    call_arguments[k] = call_arguments[k].trim();
                }
                for (k = 0; k < fd.argumentNames.length; k++) {
                    body = body.replaceAll(fd.argumentNames[k], call_arguments[k]);
                }
                expansion = expansion.replaceAll(matches[j], body);
            }
        }
        if (!done) {
            return this._expandFormula(expansion, functionDefinitions, numRecursion + 1, maxRecursion);
        }
        return expansion;
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