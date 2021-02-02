// Time-bound cache for long running processes.
// Note that a Node process will wait for the timeouts to fire
// when exiting. Use clear() to avoid hanging the process.

module.exports = {
    cache: (mk, ttl = 5000) => {
        const table = {};

        const release = (key) => {
            if (table[key]) {
                clearTimeout(table[key].alarm);
                delete table[key];
            }
        };

        const clear = () => {
            Object.keys(table).forEach(release);
        };
        
        const get = (key) => {
            if (!table[key]) {
                table[key] = {
                    value: mk(key),
                };
            }

            // Restart timer on each access.
            clearTimeout(table[key].alarm);
            table[key].alarm = setTimeout(() => clear(key), ttl);

            return table[key].value;
        };    
        
        return {
            get,
            release,
            clear,
        };
    },
};

