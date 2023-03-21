function Counter(array) {
    var count = new Map();
    array.forEach(val => count.set(val, (count.get(val) || 0) + 1));
    return count;
}

function compareMaps(map1, map2) {
    let testVal;
    if (map1.size !== map2.size) {
        return false;
    }
    for (let [key, val] of map1) {
        console.log(key, val)
        testVal = map2.get(key);
        // in cases of an undefined value, make sure the key
        // actually exists on the object so there are no false positives
        if (testVal !== val || (testVal === undefined && !map2.has(key))) {
            return false;
        }
    }
    return true;
}

/**
 * Compare two lists of strings, return true if they have the same elements
 */
function compareLists(a, b) {
    if (a.length != b.length) {
        return false;
    }
    
    var seen = {};
    a.forEach(function(v) {
        var key = (typeof v) + v;
        if (!seen[key]) {
            seen[key] = 0;
        }
        seen[key] += 1;
    });

    return b.every(function(v) {
        var key = (typeof v) + v;
        if (seen[key]) {
            seen[key] -= 1;
            return true;
        }
        // not (anymore) in the map? Wrong count, we can stop here
    });

}