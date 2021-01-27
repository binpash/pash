

```sh
cat A B | wc
# should be the same as the following
python wc.py <(cat A | wc) <(cat B | wc)
```


Aggregators are not complete, yet.

* need to take care of padding
* nee to make them cross-platform
* ideally, make them streaming-based


