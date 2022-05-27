
// These constants are picked arbitrarily, but they are big enough to explore adequate behaviors.

// The bound of the daemon input channel
#define D_IN_CHANNEL_BOUND 5
// The maximum number of fragments that will be encountered
#define NPROC 20

// To verify run: `spin -run algorithm.pml`
// Verified using `Spin Version 6.5.1 -- 3 June 2021`
// To check version, run `spin -V`


mtype = { compile, exit, done, success, fail };

// Unclear if we want bound on channel 
chan d_in = [D_IN_CHANNEL_BOUND] of { mtype, bool};
chan r_in = [D_IN_CHANNEL_BOUND] of { mtype, bool};

byte running[2];
byte running_unsafe;

active [1] proctype daemon()
{	
    // This has two cells, each indicating a different file in the file system,
    //   and if a process is currently writing to it.
    byte deps[2];
    
    bool compile_dep;
    bool dep;

    bool suc;

    xr d_in;
    xs r_in;

loop:	do
	:: d_in?done(dep) ->
        goto end
    :: d_in?exit(dep) ->
        deps[dep]--
    :: d_in?compile(compile_dep) ->
        // Non-determinism of compilation succeed
        if
        :: (true) ->
            // Compilation success 
            suc = 1;

            // Wait for my dependency
            do
            :: (deps[compile_dep] == 0) ->
                // No dependency is running
                // Increase the dependencies that are running
                deps[compile_dep] ++;
                r_in!success(compile_dep);
                break
            :: d_in?exit(dep) ->
                deps[dep] --
            od
        :: (true) -> 
            // Compilation failed
            suc = 0;

            // Wait for all
            do
            :: (deps[0] + deps[1] == 0) ->
                // Nothing is running
                r_in!fail(compile_dep);
                break
            :: d_in?exit(dep) ->
                deps[dep] --
            od
		fi
    od
end: true
}

active [1] proctype runtime()
{	
    xr r_in;
    
    // The random dep of the incoming process
    bool dep, d2;

    byte nproc;
    nproc = NPROC;

loop:
    // Non deterministic incoming fragment
    if
    :: (nproc > 0) -> dep = 0;
    :: (nproc > 0) -> dep = 1;
    :: (true) -> goto end;
    fi

    // Reduce NPROC to finish at some point
    nproc --;

    // Send a compilation request
    d_in!compile(dep);

    if 
    :: r_in?success(d2) ->
        // The response must come for the same
        assert (dep == d2);
        run Proc(dep);
    :: r_in?fail(d2) ->
        // The response must come for the same
        assert (dep == d2);

        // Run the process here
        assert (running[0] + running[1] + running_unsafe == 0);
        running_unsafe ++;
        assert (running_unsafe == 1);
        running_unsafe --;
        
        // Let the daemon know of exit
        // NOTE: This is wrong! The daemon doesn't need to know about unsafe.
        // d_in!exit(dep)
    fi

    goto loop
end: d_in!done(0)
}

proctype Proc(bool dep)
{
    running[dep] ++;
    assert (running[dep] == 1);
    running[dep] --;
    d_in!exit(dep)
}

ltl invariant {(<> daemon@end) && (<> runtime@end)}