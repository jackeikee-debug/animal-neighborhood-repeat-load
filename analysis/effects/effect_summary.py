"""Frozen descriptive paired metrics; no inferential p, pooling or CI."""
import math,statistics
def number(x):
 if x is None or x in ('NA','', 'None'):return None
 y=float(x)
 if not math.isfinite(y):raise ValueError('nonfinite effect value is not NA')
 return y
def paired(rows):
 rows=list(rows);pairs=[(r['old'],r['new']) for r in rows if r['old'] is not None and r['new'] is not None];d=[b-a for a,b in pairs];absolute=[abs(x) for x in d]
 return {'records':len(rows),'finite_pairs':len(pairs),'baseline_NA':sum(r['old'] is None for r in rows),'new_NA':sum(r['new'] is None for r in rows),'median_abs_delta_beta':statistics.median(absolute) if d else None,'max_abs_delta_beta':max(absolute) if d else None,'median_signed_delta_beta':statistics.median(d) if d else None,'min_signed_delta_beta':min(d) if d else None,'max_signed_delta_beta':max(d) if d else None,'strict_sign_changes':sum(a*b<0 for a,b in pairs),'zero_in_pair':sum(a==0 or b==0 for a,b in pairs),'inference':'descriptive conditional matched species pairs, not independent replicate tests or confidence interval'}
