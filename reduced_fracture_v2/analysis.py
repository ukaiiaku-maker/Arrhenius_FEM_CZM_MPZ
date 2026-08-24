from __future__ import annotations


class SharedEventAvalancheAnalysis:
    """Backend-neutral semantic grouping; it does not impose a time cutoff."""

    @staticmethod
    def group(transactions):
        groups=[]; current=[]
        for row in transactions:
            if current and bool(row.get("reload_before_event",False)):
                groups.append(current); current=[]
            current.append(row)
        if current: groups.append(current)
        return [{
            "physical_avalanche_index":i,
            "first_event_transaction_index":g[0]["event_transaction_index"],
            "last_event_transaction_index":g[-1]["event_transaction_index"],
            "event_count":len(g),
            "extension_m":sum(float(x["event_extension_m"]) for x in g),
            "right_censored_at_target":bool(g[-1].get("right_censored_at_target",False)),
        } for i,g in enumerate(groups)]
